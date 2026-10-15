//! Stable counting sort into reused buffers.

/// Stable sort of `src` into `dst` (same length) by `key(x)` in `0..n_keys`:
/// O(len + `n_keys`). `counts` is scratch. Equal keys keep their order in `src`.
pub(crate) fn counting_sort<T: Copy>(
    src: &[T],
    dst: &mut [T],
    counts: &mut Vec<usize>,
    n_keys: usize,
    key: impl Fn(&T) -> usize,
) {
    counts.clear();
    counts.resize(n_keys + 1, 0);
    for x in src {
        if let Some(c) = counts.get_mut(key(x) + 1) {
            *c += 1;
        }
    }
    for k in 1..counts.len() {
        counts[k] += counts[k - 1];
    }
    for x in src {
        if let Some(at) = counts.get_mut(key(x))
            && let Some(slot) = dst.get_mut(*at)
        {
            *slot = *x;
            *at += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn equals_a_stable_comparison_sort(v in proptest::collection::vec((0_usize..7, any::<u32>()), 0..200)) {
            let (mut dst, mut counts) = (v.clone(), Vec::new());
            counting_sort(&v, &mut dst, &mut counts, 7, |x| x.0);
            let mut want = v.clone();
            want.sort_by_key(|x| x.0);
            prop_assert_eq!(dst, want);
        }
    }
}
