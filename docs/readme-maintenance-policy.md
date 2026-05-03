# README maintenance policy

The README is the main scientific entry point of this repository.

## Rule

Do not reduce the README in future phases.

Updates may:

- add new sections;
- update validation status;
- add new figures;
- improve equations;
- link to expanded documentation.

Updates must not remove the following core blocks:

- thesis/project summary;
- mathematical formulation;
- MMS explanation;
- numerical verification equations;
- Case 12 status;
- source-term validation status;
- figure references;
- installation commands;
- test commands;
- examples;
- repository structure;
- Git workflow;
- roadmap.

## When the README grows too much

Move detailed derivations to `docs/`, but keep a concise explanation and link in the README.

## Curated figures

Generated figures belong in:

```text
outputs/
```

Selected README figures belong in:

```text
assets/preview/
```

Do not commit the entire `outputs/` tree.
