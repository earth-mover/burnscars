# burnscars

Rust port of the controlledburn core (`~/repos/controlledburn/cpp`).

## Upstream sources: rules for all agents and subagents

- DO NOT LOOK AT GEOS FOR INSPIRATION. Do not read it, search it, or port from it.
  GEOS is LGPL-2.1; derived code would make this crate LGPL.
- Do not use JTS either, unless the user allows it.
- For geometry algorithms, look only at exactextract (`~/repos/exactextract`, Apache-2.0).
- controlledburn (`~/repos/controlledburn`, Apache-2.0) is the code being ported; read it for the port.
- If controlledburn and exactextract disagree, exactextract wins. The oracle suite decides.
- GEOS is allowed only as a black-box test oracle (shapely in `oracle/`), never as source.
