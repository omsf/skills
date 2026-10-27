# Agent instructions

These instructions apply throughout this repository.

## Scientific and operational boundaries

- Clarify chemical identity, protonation, stereochemistry, atom ordering, units,
  and intended scientific question before adapting examples.
- Preserve the tutorial's force field and sampling choices when reproducing it.
  Record any deliberate changes, software versions, backend, and source revision.
- Legacy and prototype examples are explicitly labeled. Do not fix an API error
  by silently changing the scientific model or claiming a prototype is supported.
- Ask before paid services, large downloads, GPU runs, long simulations, or
  uploading private structures to hosted demos. Never overwrite original inputs.
- A successful example run is not evidence of convergence, predictive accuracy,
  or suitability for clinical or regulatory decisions.

## Maintain and validate

Keep one skill per directory workflow, with a matching lowercase folder/name and
a routing description. Use original instructions and small source-checked
examples, cite upstream revisions, and keep each `SKILL.md` under 500 lines. When
updating the catalog, compare the live directory with the
[README catalog](README.md#workflow-catalog) and update the corresponding skills together.

Validate frontmatter with the Agent Skills reference validator (requires `uv` and
network access on first use):

```bash
for skill_file in */SKILL.md; do
    uvx --from skills-ref agentskills validate "${skill_file%/SKILL.md}" || exit 1
done
```

Also inspect code-block syntax and source links. These checks do not run the
scientific workflows: numerical validation requires the relevant environments,
input data, optional licenses, and compute budget. Full end-to-end scientific
validation is not claimed.
