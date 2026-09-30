# PRISM compartment model — E52

Source bundle for the integrated compartmental PRISM model associated with the
SV6 epoch-52 checkpoint. The checkpoint SHA-256 is
`9fbb7df5119aba87a52d437fe443d656044e3c9c868fe13b1c39b73ac375944d`.

## Start here

- [Korean overview](README_KO.md)
- [Model card](model_card/MODEL_CARD_KO.md)
- [Checkpoint identity and provenance](model_card/CHECKPOINT_MANIFEST.json)
- Authoritative E52 configurations: [source](model_card/cfg/source_config.json)
  and [resolved](model_card/cfg/resolved_config.json).

## Code map

| Directory | Contents |
|---|---|
| `src/kmlee_bam/model/` | Encoder, common/personal paths, compartmental threshold branch and decoder |
| `src/kmlee_bam/modules/` | Module tokenizer and attention components |
| `src/kmlee_bam/objectives/` | Loss definitions |
| `src/kmlee_bam/training/` | Integrated curriculum, training loop, pathology-rank and generator selection |
| `scripts/`, `configs/`, `tests/` | Historical scripts, configurations and checks |
| `preprocessing_reference/` | Supplemental preprocessing code and environment reference |

The central compartmental implementation is
[`precision_medicine.py`](src/kmlee_bam/model/precision_medicine.py).

## Verify this copy

```bash
python3 scripts/verify_bundle.py
```

This checks file hashes, Python syntax and the E52 checkpoint identity without
loading data or starting training. With the required dependencies installed,
also check the active model/training imports:

```bash
python3 scripts/verify_bundle.py --import-check
```

## Reproduction scope

The source was copied from the E52 runtime directory on 2026-09-23. That copy
date does not establish bitwise identity with every historical training file.
The original file manifest is retained under `provenance/`.

This publication copy excludes Python caches and temporary source backups.
Model and training source files are preserved byte for byte. The old
`verify_bundle.py` checked a different rank2 E20 bundle; it is archived under
`provenance/` and replaced here by an E52 source-bundle verifier.

Raw donor data, learned checkpoints and large external Zarr/NPZ artifacts are
not included. Historical configurations retain their original absolute paths;
they require the corresponding data, reference artifacts and recovery
checkpoint before they can be used. The environment YAML in
`preprocessing_reference/` is a reference, not a historical E52 lockfile.
This is a source and review bundle, not a self-contained training reproduction.

The model card records consumed-test reuse and unresolved depth coupling.
The compartmental branch is a computational analogy to dendritic nonlinearities;
its presence does not establish that it independently improves biological recovery.
