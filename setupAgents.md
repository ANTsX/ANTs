# Setup Routed Agent Workflow

## Goal

Bootstrap a compact, repo-specific agent workflow like the one in this repository.

You are not copying another project blindly. You must inspect the target repository, infer how work is actually organized, and then create a routed instruction system that future LLM agents can follow with minimal prompt overhead.

The end state is:

- a short `AGENTS.md` entrypoint
- a `.agents/workflows/` router layer
- a `.agents/references/` layer with the smallest useful semantic docs
- routing rules that tell future agents what to read first based on query content
- a root `coding.md` generated from this document for reusable coding behavior
- a root `.gitignore` that excludes temporary/generated local files from version control
- coding/editing/notebook policies tied to real ownership boundaries in the repo
- lightweight handoff and change-log conventions

If the repository already has some of these files, treat this as a refactor or completion pass, not a rewrite. Preserve stable paths and names unless there is a strong reason to change them.

## Core Principles

- Derive from the repo. Do not invent architecture that the project does not actually have.
- Keep the entrypoint short. Push detail into routers and reference docs.
- Route by ownership and semantics, not by vague topic labels alone.
- Apply `coding.md` as the baseline coding behavior for every task; treat workflow docs as repo-specific overlays.
- Read the smallest relevant doc first.
- Open owning modules before opening large notebooks, scripts, or pipelines.
- Fix bugs at the narrowest owning layer.
- Preserve canonical outputs, filenames, variable names, and stage semantics unless a deliberate migration requires otherwise.
- If notebooks exist, keep notebooks orchestration-thin and keep reusable logic in package/module owners.
- Treat wrapper scripts and CLI entrypoints as wrapper behavior owners, not business-logic authority.
- Treat data-table, export, cache, and figure semantics as owned by the stage that writes them, not by downstream consumers.
- Require validation after edits. Do not claim success from static reasoning alone.
- Use change logs for real handoff value, not for noise.
- Optimize future prompts for brevity. The docs should let a user say only the task, target, and special constraint.

## Phase 1: Inspect The Repository

Before writing any workflow docs, inspect the repo to determine its real structure.

At minimum:

1. Read top-level docs and manifests.
   - `README*`
   - `pyproject.toml`, `package.json`, `Cargo.toml`, or equivalent
   - test config
   - existing architecture notes
2. Inventory top-level execution surfaces.
   - notebooks
   - scripts
   - CLIs
   - services/apps
   - pipelines/jobs
3. Inventory ownership layers.
   - reusable package/library code
   - orchestration layers
   - plotting/rendering layers
   - wrapper/tooling layers
   - tests and smoke checks
4. Identify canonical outputs and where they are written.
   - tables
   - caches
   - figures
   - reports
   - intermediate artifacts
5. Identify public functions or stable entrypoints that notebooks/scripts call.
6. Identify the smallest validation surface for each major workflow.
7. Identify any mixed migration state or legacy patterns that still exist but should not be treated as authoritative.

Do not start by opening large notebooks or long scripts end to end. Start from repo structure, then owner modules, then only the specific large regions you need.

## Phase 2: Define Workflow Profiles

Create a top-level router plus one or more workflow profiles.

A workflow profile is a group of tasks that share:

- the same primary entrypoint(s)
- similar stage semantics
- similar ownership modules
- similar validation surface
- similar handoff/change-log needs

Use as few profiles as possible.

Create a new profile only when the workflow genuinely needs its own:

- stage map
- routing table
- semantic references
- validation cadence
- recent-changes log and remaining-work tracker

If two notebooks or pipelines share the same owners and stage semantics, keep them in one profile.

## Phase 3: Write `AGENTS.md`

Create a short root `AGENTS.md` that future agents will read first.

It should contain:

1. A one-paragraph explanation that the repo uses a routed instruction system under `.agents/`.
2. A required startup order.
3. Repo-wide non-negotiable rules.
4. A reference-file list with one-line descriptions.
5. A scope note explaining that details live under `.agents/`.
6. A note that `coding.md` is always loaded as baseline behavior before workflow-specific routing.

Use this shape:

- `# Agent Entrypoint`
- `## Required startup order`
- `## Non-negotiable repo rules`
- `## Reference files`
- `## Scope note`

The startup order should usually be:

1. Open `coding.md` for baseline coding behavior.
2. Open the top-level router.
3. Follow its dispatch table.
4. Read the smallest relevant reference doc.
5. Open symbol/stage maps only if needed.
6. Open owning modules before large notebooks/scripts.
7. Open large notebook/script regions only when owner-module context is insufficient.

## Phase 4: Write Routers

Create:

- one thin top-level router in `.agents/workflows/`
- one router per workflow profile

### Top-Level Router

The top-level router should be a dispatcher only.

Use this shape:

- title: `# <Project> Router`
- `Purpose`
- `Use this file when`
- `## Read this first`
- `## Workflow profile dispatch`
- `## Cross-workflow invariants`
- `## Compact scaling rule`

Its dispatch table should map primary targets to profile routers.

Its cross-workflow invariants should include only truly repo-wide rules, such as:

- prefer package/module edits over orchestration-layer edits
- preserve canonical outputs and stage semantics
- do not use wrappers as business-logic authority
- fix semantics at the writer stage, not a downstream consumer
- apply `coding.md` baseline behavior, then layer repo-specific rules from routers/references

Its scaling rule should explain when to add an entrypoint to an existing profile versus when to create a new profile.

### Profile Routers

Each profile router should be compact and operational.

Use this shape:

- title: `# <Profile> Router`
- `Purpose`
- `Use this file when`
- `## Read order`
- `## Task routing table`
- `## Ownership guidance`

Each profile router should:

- define the read order for that profile
- map common query types to the smallest relevant reference docs
- name the owning module/package/script layers
- name the correct workflow-specific recent-changes log and remaining-work tracker
- treat `coding.md` as global behavior and avoid duplicating generic coding guidance in profile docs

The task routing table is the core value. Route by likely semantic owner.

Examples of task categories:

- stage ownership or setup
- registration/transforms/orientation
- canonical outputs or schema confusion
- figures and plotting
- cache/rerun behavior
- response/classification semantics
- notebook refactor ownership
- CLI wrapper behavior
- migration continuation or resume point

Adapt those categories to the target repo. Do not force notebook-specific terms into a repo that does not use notebooks.

## Phase 5: Write Reference Docs

Create only the smallest set of reference docs that future agents will repeatedly need. Prefer several short docs over one giant catch-all doc.

Common reference docs:

- `stage-map.md` or `<profile>-stage-map.md`
- `symbol-index.md`
- `current-state.md`
- `refactor-rules.md`
- `refactor-loop-policy.md`
- `recent-changes.md`
- `recent-changes-<profile>.md`
- `remaining-work.md`
- `remaining-work-<profile>.md`

Create additional semantic docs only when they correspond to real repo concerns, such as:

- `canonical-outputs.md`
- `data-contracts.md`
- `activity-semantics.md`
- `figure-rules.md`
- `cache-rerun-policy.md`
- `scientific-policy.md`
- `api-semantics.md`
- `artifact-ownership.md`

Do not create empty placeholders. If a category is not meaningful in the target repo, omit it.

### Stage Map

Use a stage map when the workflow has a stable ordered flow through a notebook, script, or pipeline.

Include:

- purpose
- when to use the file
- end-to-end stages in order
- key outputs by stage
- concept ownership
- navigation notes pointing to symbol index or semantic docs

If the repo uses tagged notebook cells, record the cell tags. If it uses scripts or pipeline steps instead, record function names, targets, commands, or modules.

### Symbol Index

Create a compact symbol index for the extracted public surface the notebooks/scripts should call.

Include:

- trusted implementation patterns
- module-by-module public symbols
- one-line ownership notes for each public function/class

Focus on public, notebook-callable, script-callable, or stage-owning functions. Do not try to document every private helper.

### Current State

Use `current-state.md` to record mixed migration state, known semantic caveats, or temporary practical warnings.

Put unstable migration notes here instead of spreading them across all routers.

### Refactor Rules

Use `refactor-rules.md` for ownership and edit-scope rules.

Include rules such as:

- where orchestration belongs
- where reusable logic belongs
- where figure construction belongs
- where wrapper behavior belongs
- no new local helpers in orchestration layers when reusable code should be extracted
- fix at the narrowest owner
- do not patch downstream consumers to compensate for upstream semantic bugs
- preserve canonical outputs and stage order
- prefer package-level smoke checks over ad hoc assertions

Make these rules repo-specific. Replace notebook language with script/pipeline language if that better fits the project.

### Refactor Loop Policy

Use `refactor-loop-policy.md` to keep agents moving through an ownership slice instead of stopping after one superficial cleanup.

Define:

- the default working unit
- the required sub-pass loop
- keep-going rules
- valid stop conditions
- invalid stop conditions
- handoff requirements when stopping mid-slice
- validation expectations

The key idea is: one cleaned-up function or one cleaned-up cell is not automatically a valid stopping point if the remaining work is clearly inside the same owner and validation surface.

### Recent Changes Logs And Remaining Work Trackers

Create:

- one recent-changes index file that routes to workflow-specific change logs
- one rolling recent-changes log per workflow profile
- one remaining-work index file that routes to workflow-specific remaining-work trackers
- one rolling remaining-work tracker per workflow profile

The recent-changes index file should say which workflow change log to open. The remaining-work index file should say which workflow tracker to open.

Each workflow recent-changes log should include an append-only template for completed meaningful work, such as:

- date and short label
- slice goal
- passes completed in this session
- what changed
- rerun implications
- validation performed

Log every meaningful change, not only refactors or unfinished work. A meaningful change is one that affects behavior, public workflow semantics, canonical outputs, validation expectations, ownership boundaries, or future rerun/debugging decisions.

Use these logs for real handoff value, not for noise. Do not log purely mechanical formatting, typo-only edits, or regenerated artifacts unless they change semantics or future agent behavior.

Keep unresolved breakage and follow-up work out of recent-changes logs. Each workflow remaining-work tracker should include an append-only or status-updated template such as:

- date and short label
- current status
- what remains broken
- remaining in-slice work
- next likely breakpoint
- blocking context
- validation or rerun needed

Use remaining-work trackers for actionable unresolved work, not for historical narration. When a tracked item is resolved, mark it resolved in the tracker and record the completed fix in the workflow recent-changes log.

## Phase 6: Derive Repo Rules From Real Ownership

The workflow is only useful if its policies reflect the real codebase.

Infer and encode:

- which layer owns reusable logic
- which layer owns orchestration
- which layer owns rendering/figures
- which layer owns wrappers or CLI glue
- which stage writes authoritative outputs
- which downstream consumers should never redefine upstream semantics
- which legacy patterns still exist but are not authoritative

If the repo contains notebooks, use notebook-specific rules like:

- notebooks remain orchestration-thin
- reusable logic lives in modules/packages
- no new notebook-local helpers when extraction is appropriate
- preserve notebook stage order, cell tags, outputs, and legacy variable names unless intentionally migrating them

If the repo does not contain notebooks, apply the same principle to:

- orchestration scripts
- pipeline entrypoints
- task runners
- service handlers

The policy should describe the actual ownership boundary, not a Python-specific ideology.

## Phase 7: Define Future Read Order

Your docs should make future task handling predictable.

The default future read order should usually be:

1. `AGENTS.md`
2. `coding.md`
3. top-level router
4. profile router
5. smallest relevant semantic reference
6. stage map or symbol index only if still needed
7. owning module/package/script
8. large notebook/script region only if owner context is insufficient
9. wrappers only for wrapper behavior

Write this ordering into the routers and entrypoint explicitly.

## Phase 8: Define Task Routing By Query Content

Future agents should classify each task along these axes:

- primary target: which notebook, script, CLI, service, or module
- workflow profile: which stage family or pipeline
- semantic owner: outputs, matching, plotting, caching, API, orchestration, etc.
- ownership layer: package/module, orchestration, wrapper, tests
- validation surface: smoke test, unit test, contract check, rerun stage, downstream consumer

Then the agent should open the smallest relevant reference doc first.

This means your task-routing tables should map query content to reference docs, not just to source files.

Bad route:

- "figure bug -> open notebook"

Better route:

- "figure bug -> open figure rules, then stage map if needed, then figure owner module, then notebook/script caller only if necessary"

## Phase 9: Validation And Maintenance Rules

Encode the following operating rules into the workflow:

- run the smallest relevant smoke, contract, or regression check after edits
- when a writer stage or canonical output changes, verify the first downstream consumer
- when changes affect plots, figures, rendered reports, notebooks, GUIs, dashboards, or frontend views, verify the rendered output with screenshots or artifact inspection; check for blank renders, wrong data, clipping, overlap, unreadable labels, broken layout, and incorrect visual state
- for web UIs, prefer automated browser and screenshot tooling such as Playwright when available
- do not declare success from static reasoning alone
- if public behavior changes, update the relevant reference doc
- if public function surface changes, update `symbol-index.md`
- if stage/output ownership changes, update the stage map and `current-state.md`
- append the workflow-specific recent-changes log for every meaningful change, including completed work
- if meaningful work stops with remaining breakage or follow-up work, update the workflow-specific remaining-work tracker instead of burying unresolved work in the change log
- when tracked remaining work is resolved, mark it resolved in the remaining-work tracker and log the completed fix in recent-changes
- if baseline coding guidance changes, regenerate `coding.md` and keep `AGENTS.md` startup order aligned

If the repo lacks tests, define the smallest practical validation surface anyway:

- a targeted script run
- a notebook stage rerun
- a CLI dry run
- a golden-file comparison
- a static contract check
- a screenshot or rendered-artifact check for visual outputs

## Phase 9a: Produce `coding.md`

Create a root `coding.md` alongside `AGENTS.md`.

Purpose:

- keep generic coding behavior separate from repo-specific routing docs
- give future prompts a stable short file to include when they want coding behavior only
- keep one canonical source of truth for those guidelines

Rules:

- `setupAgents.md` is the source of truth for `coding.md`
- store the exact `coding.md` body in the appendix block between the `coding.md:start` and `coding.md:end` markers
- generate `coding.md` from that block instead of maintaining two hand-edited copies
- if you change the coding guidance, regenerate `coding.md` in the same pass
- keep `coding.md` generic unless a rule truly belongs in every repo
- keep repo-specific ownership and routing rules in `AGENTS.md` and `.agents/`, not in `coding.md`
- require `AGENTS.md` and the top-level router to point to `coding.md` as baseline behavior

The file should stay short and operational. It is meant to capture stable coding behavior such as:

- think before coding
- simplicity-first bias
- surgical changes
- goal-driven execution

## Phase 9b: Produce Or Update `.gitignore`

Create or update a root `.gitignore` so generated local artifacts are not tracked.

At minimum, include ignore rules for common temporary files relevant to the repo, such as:

- Python bytecode and caches: `__pycache__/`, `*.py[cod]`, `.pytest_cache/`, `.mypy_cache/`, `.ruff_cache/`
- virtual environments: `.venv/`, `venv/`, `env/`
- build and packaging outputs: `build/`, `dist/`, `*.egg-info/`
- editor and OS noise: `.DS_Store`, `.idea/`, `.vscode/`
- runtime scratch files that the repo produces locally, such as logs, temporary outputs, or cache directories

Keep `.gitignore` repo-specific. Do not ignore canonical outputs, checked-in fixtures, or generated artifacts that the project intentionally versions.

## Phase 10: Quality Bar For The Docs

The docs should be:

- short
- operational
- specific
- stable
- non-redundant

Avoid:

- giant architecture essays
- duplicated rules in every file
- vague "best practices" with no ownership mapping
- docs that require reading the whole repo before being useful
- empty placeholders with only TODO text

When uncertain:

- prefer a narrow, defensible statement
- record open ambiguity in `current-state.md`
- avoid pretending that guessed semantics are settled

## Deliverables Checklist

Before you consider the setup complete, make sure you have produced:

- `AGENTS.md`
- `coding.md`
- `.gitignore` with temporary, cache, build, editor, and OS-noise patterns appropriate for the repo
- `.agents/workflows/<top-level-router>.md`
- one profile router per real workflow profile
- stage maps for workflows that need ordered-stage navigation
- `symbol-index.md`
- `refactor-rules.md`
- `refactor-loop-policy.md`
- `current-state.md` if there is active migration or semantic ambiguity
- `recent-changes.md`
- workflow-specific `recent-changes-*.md` files
- `remaining-work.md`
- workflow-specific `remaining-work-*.md` files
- the smallest semantic reference docs needed for repeated routing decisions
- explicit references in `AGENTS.md` and top-level router that `coding.md` is required baseline guidance

Then sanity-check the workflow by routing a few representative task prompts. For each prompt, verify that the docs make it obvious:

1. which router to open
2. which reference doc to read first
3. which owner layer to edit
4. which validation to run
5. which log/doc to update if the change affects public behavior
6. how `coding.md` baseline behavior applies without duplicating generic guidance across repo-specific docs

## Final Instruction

Do not optimize for a perfect doc tree on the first pass. Optimize for a compact routed system that future agents can actually use.

A good setup makes future prompts smaller, owner boundaries clearer, notebook/script edits rarer, and validation expectations explicit.

## Appendix: `coding.md` Source

Edit this block, then regenerate `coding.md`.

<!-- coding.md:start -->
> Generated from `setupAgents.md` by `scripts/sync_coding_doc.py`. Do not edit directly.

Behavioral guidelines to reduce common LLM coding mistakes. Merge with project-specific instructions as needed.

Tradeoff: These guidelines bias toward caution over speed. For trivial tasks, use judgment.

1. Think Before Coding

Don't assume. Don't hide confusion. Surface tradeoffs.

Before implementing:

State your assumptions explicitly. If uncertain, ask.
If multiple interpretations exist, present them - don't pick silently.
If a simpler approach exists, say so. Push back when warranted.
If something is unclear, stop. Name what's confusing. Ask.

2. Simplicity First

Minimum code that solves the problem. Nothing speculative.

No features beyond what was asked.
No abstractions for single-use code.
No "flexibility" or "configurability" that wasn't requested.
No error handling for impossible scenarios.
If you write 200 lines and it could be 50, rewrite it.
Ask yourself: "Would a senior engineer say this is overcomplicated?" If yes, simplify.

3. Surgical Changes

Touch only what you must. Clean up only your own mess.

When editing existing code:

Don't "improve" adjacent code, comments, or formatting.
Don't refactor things that aren't broken.
Match existing style, even if you'd do it differently.
If you notice unrelated dead code, mention it - don't delete it.

When your changes create orphans:

Remove imports/variables/functions that YOUR changes made unused.
Don't remove pre-existing dead code unless asked.
The test: Every changed line should trace directly to the user's request.

4. Goal-Driven Execution

Define success criteria. Loop until verified.

Transform tasks into verifiable goals:

"Add validation" -> "Write tests for invalid inputs, then make them pass"
"Fix the bug" -> "Write a test that reproduces it, then make it pass"
"Refactor X" -> "Ensure tests pass before and after"
"Fix plot/UI layout" -> "Capture the rendered output and verify the screenshot shows correct data, labels, spacing, and state"

For multi-step tasks, state a brief plan:

1. [Step] -> verify: [check]
2. [Step] -> verify: [check]
3. [Step] -> verify: [check]

Strong success criteria let you loop independently. Weak criteria ("make it work") require constant clarification.
<!-- coding.md:end -->
