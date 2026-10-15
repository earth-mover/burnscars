# Task runner. `just` lists the recipes.

default:
    @just --list

build:
    cargo build --all-targets

check:
    cargo check --all-targets --all-features

lint:
    cargo clippy --all-targets --all-features -- -D warnings

format *args:
    cargo fmt --all {{ args }}

alias fmt := format

format-toml *args:
    taplo fmt {{ args }}

# All tests (nextest). Doctests run in `doctest`.
test *args:
    cargo nextest run --no-tests=pass {{ args }}

# In-source unit tests only, filtered (e.g. `just test-lib grid::`).
test-lib *args:
    cargo nextest run --lib --no-tests=pass {{ args }}

# One integration test file, filtered (e.g. `just test-one invariants offset`).
test-one name *args:
    cargo nextest run --test {{ name }} --no-tests=pass {{ args }}

doctest:
    cargo test --doc

bench *args:
    cargo bench --bench burn {{ args }}

# Regenerate the oracle fixtures (needs uv).
oracle *args:
    uv run --project oracle python oracle/generate.py --out fixtures/oracle {{ args }}

# Differential fuzz against the live oracle (needs uv).
fuzz-oracle:
    cargo nextest run --test oracle --run-ignored only

check-deps:
    cargo deny check

ci:
    just format --check
    just lint
    just test
    just doctest

check-all: ci
    just check-deps
    just format-toml --check
