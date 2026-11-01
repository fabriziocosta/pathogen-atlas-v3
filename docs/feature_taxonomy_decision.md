# Feature taxonomy decision

## Decision

Organize the V3 feature catalogue as a three-level hierarchy:

1. **Dimension** — the broad kind of information being measured.
2. **Category** — a coherent group of related concepts within that dimension.
3. **Feature** — one versioned question assessed against the title and abstract.

The hierarchy lives in `configs/features.yaml`, so it can be inspected and edited alongside each feature's ID, definition, applicable organism groups, question version, evidence instructions, criteria, and review status. Every feature uses the shared evidence frame:

> In this abstract, is this specific claim explicitly reported about the named species in a human disease context: “[feature claim]”?

Each feature also has an explicit `question_claim`, plus `inclusion_examples` and `exclusion_examples`. The examples are schematic claim-specific patterns, not quotations or claims about particular publications. They clarify what counts as evidence for that feature and when similar wording is insufficient. They are included in Jev's instructions. `configs/feature_profiles.yaml` remains the explicit ordered list of features to score for each organism group. The Python loader validates and flattens the hierarchy for scoring; the Explorer consumes the same hierarchy to group feature summaries and show the question and examples.

## Dimensions

| Dimension | What belongs here |
| --- | --- |
| **Pathogen biology** | Organism traits, infection course, and response to antimicrobial treatment. |
| **Transmission** | Routes and patterns by which the organism or infection spreads between hosts. |
| **Reservoirs and exposure** | Animal or environmental reservoirs and settings associated with exposure. |
| **Clinical presentation** | Human disease manifestations, affected systems, diagnostic findings, and co-occurring conditions. |
| **Host factors** | Patient conditions, behaviours, immune status, and susceptibility. |
| **Healthcare context** | Care settings, devices, procedures, and clinical products associated with infection. |
| **Epidemiology and public health** | Population patterns, outbreak occurrence, public-health classification, prevention, and control. |
| **Abstract reporting** | Study or case descriptors reported in an abstract, such as location, age, occupation, or activity. These describe the report's content, not a property of the organism. |

Current category names and membership are visible directly under each dimension in the YAML. The catalogue is a working structure, not a claim that every inherited concept has already passed expert review.

## Placement and editing rules

- Assign each feature one primary dimension and one category. Avoid duplicate questions as a way to represent cross-cutting concepts; if needed, add searchable tags later without duplicating the feature.
- Classify by what the question asks the abstract to establish, not by which team owns it or which organism group uses it.
- Keep a feature specific enough to answer consistently from an abstract. Separate questions that combine distinct routes, hosts, syndromes, or mechanisms.
- Use explicit human-disease and named-organism context where relevant. A positive answer means the abstract supports the feature under its wording; it does not establish causality beyond the source evidence.
- Distinguish “not mentioned or unsupported” from “not applicable.” Profiles determine applicability; an omitted feature must not be interpreted as a negative answer.
- Move a feature between dimensions or categories without changing its ID or question version when its meaning is unchanged. If the wording, criteria, scope, or interpretation changes, issue a new question version so cached decisions are rescored.
- Keep questionable, overlapping, or poorly defined inherited features visible and marked for review. Do not silently turn a broad or ambiguous question into a confident biological claim.

## Review implications

Some current questions still need subject-matter review for overlap or scope. Examples include animal group features versus zoonotic transmission, airborne as a route versus air as an environmental association, overlapping immune-status questions, and disease or diagnostic concepts grouped near clinical presentation. “Priority species” also depends on an external authority and date. These remain inspectable so they can be clarified, split, merged, or retired deliberately. Any semantic edits must update the question version and prompt/evaluation records before broad rescoring.

## Versioning

`schema_version: 2` identifies the nested YAML structure. The catalogue is now `multi-domain-features-v3.2`; this changes the feature prompt because the shared evidence frame and examples are now part of Jev's instructions. Previous v3.1 feature decisions therefore do not count as cached v3.2 decisions and will be rescored on the next Phase 3 run. Domain profiles continue to control the applicable feature set and its order. The Phase 3 cost estimate was raised to account for the larger feature instructions.

## V2 bacteria initialization

V3's 73-feature bacteria profile has the same feature IDs and definitions as V2's 73-question set. To avoid repeating this work, V3 may initialize its bacteria decision ledger with V2 pathogen-gate and feature probabilities, using the V3 run IDs and question/profile versions so matching decisions act as V3 cache entries. This is an explicit compatibility assumption: V2 used earlier prompt instructions, while V3 added the shared evidence frame and inclusion/exclusion examples. Each imported row therefore retains `decision_source: v2_compatibility_seed` and the original V2 question version, run ID, and input hash; feature rows also retain their V2 feature-set version. V2 did not record a separate gate profile version. They are not represented as newly returned V3 provider responses.

Only decisions from the pinned `typesafe/jev-1.14` model are eligible. V2 and V3 must have the same selected species–PMID pair, matched aliases, title, and abstract; changed or absent inputs remain pending for V3 Jev scoring. The import covers the V2 disease gate as well as feature probabilities, so the matching V2 confidence is adopted for the V3 gate threshold. Disable `migration.seed_v2_bacteria_decisions_as_v3` in `configs/default.yaml` to require fresh V3 decisions instead.
