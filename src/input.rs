//! Input normalization: geo-traits geometry into a reused buffer.
//!
//! Rules from the spec "Input normalization". Orientation by controlledburn `signed_area`
//! (Apache-2.0, Copyright 2025 Michael Sumner), commit 1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558,
//! `cpp/src/scanline.cpp`.

use std::ops::Range;

use geo_traits::{
    CoordTrait as _, GeometryCollectionTrait as _, GeometryTrait, GeometryType,
    LineStringTrait as _, LineTrait as _, MultiLineStringTrait as _, MultiPointTrait as _,
    MultiPolygonTrait as _, PointTrait as _, PolygonTrait, RectTrait as _, TriangleTrait as _,
};

use crate::error::{At, BurnError, GeomKind};
use crate::geom::{coord::Coord, measures::signed_area_ccw};
use crate::grid::GeomId;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum RingRole {
    Exterior,
    Hole,
}

impl RingRole {
    pub(crate) fn sign_f32(self) -> f32 {
        match self {
            Self::Exterior => 1.0,
            Self::Hole => -1.0,
        }
    }
    pub(crate) fn sign_i32(self) -> i32 {
        match self {
            Self::Exterior => 1,
            Self::Hole => -1,
        }
    }
}

/// Closed, >= 4 coordinates, all finite, CCW. Only `InputBuf` makes one.
#[derive(Clone, Copy, Debug)]
pub(crate) struct NormRing<'s> {
    pub(crate) coords: &'s [Coord],
    pub(crate) role: RingRole,
    pub(crate) part: usize,
    pub(crate) ring: usize,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Envelope {
    pub(crate) xmin: f64,
    pub(crate) ymin: f64,
    pub(crate) xmax: f64,
    pub(crate) ymax: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Loaded {
    Empty,
    Polygons,
    Lines,
    Points,
}

#[derive(Clone, Debug)]
struct RingMeta {
    range: Range<usize>,
    role: RingRole,
    part: usize,
    ring: usize,
}

#[derive(Clone, Debug)]
struct LineMeta {
    range: Range<usize>,
    part: usize,
}

/// Reused normalized geometry. One geometry at a time.
#[derive(Debug, Default)]
pub(crate) struct InputBuf {
    coords: Vec<Coord>,
    rings: Vec<RingMeta>,
    lines: Vec<LineMeta>,
    points: Vec<Coord>,
}

impl InputBuf {
    pub(crate) fn rings(&self) -> impl Iterator<Item = NormRing<'_>> {
        self.rings.iter().filter_map(|m| {
            Some(NormRing {
                coords: self.coords.get(m.range.clone())?,
                role: m.role,
                part: m.part,
                ring: m.ring,
            })
        })
    }

    pub(crate) fn lines(&self) -> impl Iterator<Item = (usize, &[Coord])> {
        self.lines.iter().filter_map(|m| Some((m.part, self.coords.get(m.range.clone())?)))
    }

    pub(crate) fn points(&self) -> &[Coord] {
        &self.points
    }
}

impl InputBuf {
    fn clear(&mut self) {
        self.coords.clear();
        self.rings.clear();
        self.lines.clear();
        self.points.clear();
    }

    /// Validates and normalizes one geometry before it writes any record.
    pub(crate) fn load<G: GeometryTrait<T = f64>>(
        &mut self,
        g: &G,
        id: GeomId,
    ) -> Result<Loaded, BurnError> {
        self.clear();
        let at = |part: usize| At { id, part, ring: 0, vertex: 0 };
        match g.as_type() {
            GeometryType::Point(p) => self.push_point(p.coord().map(|c| c.x_y()), at(0))?,
            GeometryType::MultiPoint(mp) => {
                for (i, p) in mp.points().enumerate() {
                    self.push_point(p.coord().map(|c| c.x_y()), at(i))?;
                }
            }
            GeometryType::LineString(l) => self.push_line(l.coords().map(|c| c.x_y()), at(0))?,
            GeometryType::MultiLineString(m) => {
                for (i, l) in m.line_strings().enumerate() {
                    self.push_line(l.coords().map(|c| c.x_y()), at(i))?;
                }
            }
            GeometryType::Line(l) => {
                self.push_line(l.coords().into_iter().map(|c| c.x_y()), at(0))?;
            }
            GeometryType::Polygon(p) => self.push_polygon(p, at(0))?,
            GeometryType::MultiPolygon(m) => {
                for (i, p) in m.polygons().enumerate() {
                    self.push_polygon(&p, at(i))?;
                }
            }
            GeometryType::Rect(r) => {
                let ((x0, y0), (x1, y1)) = (r.min().x_y(), r.max().x_y());
                let ring = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)];
                self.push_ring(ring.into_iter(), at(0), RingRole::Exterior)?;
            }
            GeometryType::Triangle(t) => {
                self.push_ring(t.coords().into_iter().map(|c| c.x_y()), at(0), RingRole::Exterior)?;
            }
            GeometryType::GeometryCollection(_) => {
                return Err(BurnError::Unsupported { id, kind: GeomKind::GeometryCollection });
            }
        }
        Ok(if !self.rings.is_empty() {
            Loaded::Polygons
        } else if !self.lines.is_empty() {
            Loaded::Lines
        } else if !self.points.is_empty() {
            Loaded::Points
        } else {
            Loaded::Empty
        })
    }

    fn push_polygon<P: PolygonTrait<T = f64>>(&mut self, p: &P, at: At) -> Result<(), BurnError> {
        let Some(ext) = p.exterior() else {
            return Ok(());
        };
        let before = self.rings.len();
        self.push_ring(ext.coords().map(|c| c.x_y()), at, RingRole::Exterior)?;
        if self.rings.len() == before {
            return Ok(()); // empty exterior: skip the holes too
        }
        for (i, h) in p.interiors().enumerate() {
            self.push_ring(h.coords().map(|c| c.x_y()), At { ring: i + 1, ..at }, RingRole::Hole)?;
        }
        Ok(())
    }

    /// Copies finite coordinates; returns the range, or `None` for an empty input.
    fn push_coords(
        &mut self,
        coords: impl Iterator<Item = (f64, f64)>,
        at: At,
    ) -> Result<Option<Range<usize>>, BurnError> {
        let start = self.coords.len();
        for (vertex, (x, y)) in coords.enumerate() {
            if !(x.is_finite() && y.is_finite()) {
                self.coords.truncate(start);
                return Err(BurnError::NonFinite { at: At { vertex, ..at }, x, y });
            }
            self.coords.push(Coord::new(x, y));
        }
        Ok((self.coords.len() > start).then_some(start..self.coords.len()))
    }

    fn push_ring(
        &mut self,
        coords: impl Iterator<Item = (f64, f64)>,
        at: At,
        role: RingRole,
    ) -> Result<(), BurnError> {
        let Some(range) = self.push_coords(coords, at)? else {
            return Ok(());
        };
        let (first, last) = (self.coords[range.start], self.coords[range.end - 1]);
        if first != last {
            self.coords.push(first); // auto-close: geo-traits does not guarantee closure
        }
        let range = range.start..self.coords.len();
        if range.len() < 4 {
            // controlledburn skips such a ring silently; burnscars reports it.
            return Err(BurnError::ShortRing { at, n: range.len() });
        }
        let ring = &mut self.coords[range.clone()];
        // Orientation as controlledburn `signed_area`: https://github.com/hypertidy/controlledburn/blob/1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558/cpp/src/scanline.cpp#L42
        let area = signed_area_ccw(ring);
        if area.is_nan() || area <= 0.0 {
            ring.reverse();
        }
        self.rings.push(RingMeta { range, role, part: at.part, ring: at.ring });
        Ok(())
    }

    fn push_line(
        &mut self,
        coords: impl Iterator<Item = (f64, f64)>,
        at: At,
    ) -> Result<(), BurnError> {
        let Some(range) = self.push_coords(coords, at)? else {
            return Ok(());
        };
        if range.len() == 1 {
            return Err(BurnError::ShortLine { at });
        }
        self.lines.push(LineMeta { range, part: at.part });
        Ok(())
    }

    fn push_point(&mut self, xy: Option<(f64, f64)>, at: At) -> Result<(), BurnError> {
        let Some((x, y)) = xy else {
            return Ok(()); // POINT EMPTY
        };
        if !(x.is_finite() && y.is_finite()) {
            return Err(BurnError::NonFinite { at, x, y });
        }
        self.points.push(Coord::new(x, y));
        Ok(())
    }

    pub(crate) fn envelope(&self) -> Option<Envelope> {
        self.coords.iter().chain(&self.points).fold(None, |env, c| {
            Some(match env {
                None => Envelope { xmin: c.x, ymin: c.y, xmax: c.x, ymax: c.y },
                Some(e) => Envelope {
                    xmin: e.xmin.min(c.x),
                    ymin: e.ymin.min(c.y),
                    xmax: e.xmax.max(c.x),
                    ymax: e.ymax.max(c.y),
                },
            })
        })
    }
}

pub(crate) fn for_each_xy<G: GeometryTrait<T = f64>>(g: &G, f: &mut dyn FnMut(f64, f64)) {
    let mut each = |(x, y): (f64, f64)| f(x, y);
    match g.as_type() {
        GeometryType::Point(p) => p.coord().map(|c| c.x_y()).into_iter().for_each(&mut each),
        GeometryType::MultiPoint(mp) => {
            mp.points().filter_map(|p| p.coord().map(|c| c.x_y())).for_each(&mut each);
        }
        GeometryType::LineString(l) => l.coords().map(|c| c.x_y()).for_each(&mut each),
        GeometryType::MultiLineString(m) => {
            m.line_strings().for_each(|l| l.coords().map(|c| c.x_y()).for_each(&mut each));
        }
        GeometryType::Line(l) => l.coords().into_iter().map(|c| c.x_y()).for_each(&mut each),
        GeometryType::Polygon(p) => polygon_xy(p, &mut each),
        GeometryType::MultiPolygon(m) => m.polygons().for_each(|p| polygon_xy(&p, &mut each)),
        GeometryType::Rect(r) => [r.min().x_y(), r.max().x_y()].into_iter().for_each(&mut each),
        GeometryType::Triangle(t) => t.coords().into_iter().map(|c| c.x_y()).for_each(&mut each),
        GeometryType::GeometryCollection(gc) => gc.geometries().for_each(|g| for_each_xy(&g, f)),
    }
}

fn polygon_xy<P: PolygonTrait<T = f64>>(p: &P, each: &mut impl FnMut((f64, f64))) {
    p.exterior()
        .into_iter()
        .chain(p.interiors())
        .for_each(|r| r.coords().map(|c| c.x_y()).for_each(&mut *each));
}

#[cfg(test)]
mod tests {
    use super::*;
    use geo_types::{
        Geometry, GeometryCollection, LineString, MultiPolygon, Point, Polygon, Rect, coord,
    };

    fn ls(pts: &[(f64, f64)]) -> LineString<f64> {
        pts.iter().map(|&(x, y)| coord! { x: x, y: y }).collect()
    }

    fn load<G: GeometryTrait<T = f64>>(g: &G) -> (InputBuf, Result<Loaded, BurnError>) {
        let mut b = InputBuf::default();
        let r = b.load(g, GeomId(3));
        (b, r)
    }

    #[test]
    fn clockwise_ring_is_reversed_and_hole_has_role() {
        let p = Polygon::new(
            ls(&[(0.0, 0.0), (0.0, 4.0), (4.0, 4.0), (4.0, 0.0), (0.0, 0.0)]), // CW
            vec![ls(&[(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0), (1.0, 1.0)])],
        );
        let (b, r) = load(&p);
        assert_eq!(r.expect("ok"), Loaded::Polygons);
        let rings: Vec<NormRing<'_>> = b.rings().collect();
        assert_eq!(rings.len(), 2);
        assert!(signed_area_ccw(rings[0].coords) > 0.0);
        assert!(signed_area_ccw(rings[1].coords) > 0.0);
        assert_eq!((rings[0].role, rings[1].role), (RingRole::Exterior, RingRole::Hole));
        assert_eq!(rings[1].ring, 1);
    }

    #[test]
    fn unclosed_ring_is_closed() {
        // A triangle from Rect-free input: geo-types closes rings, so use a Triangle.
        let t = geo_types::Triangle::new(
            coord! { x: 0.0, y: 0.0 },
            coord! { x: 1.0, y: 0.0 },
            coord! { x: 0.0, y: 1.0 },
        );
        let (b, r) = load(&t);
        assert_eq!(r.expect("ok"), Loaded::Polygons);
        let ring = b.rings().next().expect("ring");
        assert_eq!(ring.coords.len(), 4);
        assert_eq!(ring.coords[0], ring.coords[3]);
    }

    #[test]
    fn rect_becomes_ccw_ring() {
        let r = Rect::new(coord! { x: 0.0, y: 0.0 }, coord! { x: 2.0, y: 1.0 });
        let (b, res) = load(&r);
        assert_eq!(res.expect("ok"), Loaded::Polygons);
        let ring = b.rings().next().expect("ring");
        assert_eq!(ring.coords.len(), 5);
        assert_eq!(signed_area_ccw(ring.coords), 2.0);
    }

    #[test]
    fn short_ring_is_an_error() {
        let p = Polygon::new(ls(&[(0.0, 0.0), (1.0, 1.0), (0.0, 0.0)]), vec![]);
        let (_, r) = load(&p);
        assert!(matches!(r, Err(BurnError::ShortRing { n: 3, .. })), "{r:?}");
    }

    #[test]
    fn non_finite_reports_the_vertex() {
        let p = Polygon::new(
            ls(&[(0.0, 0.0), (1.0, 0.0), (f64::NAN, 1.0), (0.0, 1.0), (0.0, 0.0)]),
            vec![],
        );
        let (_, r) = load(&p);
        match r {
            Err(BurnError::NonFinite { at, .. }) => {
                assert_eq!((at.id, at.ring, at.vertex), (GeomId(3), 0, 2));
            }
            other => panic!("expected NonFinite, got {other:?}"),
        }
    }

    #[test]
    fn short_line_and_collection_are_errors() {
        let l = ls(&[(1.0, 1.0)]);
        assert!(matches!(load(&l).1, Err(BurnError::ShortLine { .. })));
        let gc: GeometryCollection<f64> =
            GeometryCollection::new_from(vec![Geometry::Point(Point::new(0.0, 0.0))]);
        assert!(matches!(
            load(&gc).1,
            Err(BurnError::Unsupported { kind: GeomKind::GeometryCollection, .. })
        ));
    }

    #[test]
    fn empty_geometries_load_as_empty() {
        let p = Polygon::new(LineString::new(vec![]), vec![]);
        assert_eq!(load(&p).1.expect("ok"), Loaded::Empty);
        let mp: MultiPolygon<f64> = MultiPolygon::new(vec![]);
        assert_eq!(load(&mp).1.expect("ok"), Loaded::Empty);
    }

    #[test]
    fn envelope_covers_all_parts() {
        let mp = MultiPolygon::new(vec![
            Polygon::new(ls(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)]), vec![]),
            Polygon::new(ls(&[(5.0, 5.0), (6.0, 5.0), (6.0, 7.0), (5.0, 5.0)]), vec![]),
        ]);
        let (b, _) = load(&mp);
        assert_eq!(b.envelope(), Some(Envelope { xmin: 0.0, ymin: 0.0, xmax: 6.0, ymax: 7.0 }));
    }

    #[test]
    fn for_each_xy_visits_collections() {
        let gc: GeometryCollection<f64> = GeometryCollection::new_from(vec![
            Geometry::Point(Point::new(1.0, 2.0)),
            Geometry::LineString(ls(&[(3.0, 4.0), (5.0, 6.0)])),
        ]);
        let mut seen = Vec::new();
        for_each_xy(&gc, &mut |x, y| seen.push((x, y)));
        assert_eq!(seen, vec![(1.0, 2.0), (3.0, 4.0), (5.0, 6.0)]);
    }
}
