# Feature taxonomy decision

## Decision

Organize the V3 feature catalogue as a three-level hierarchy:

1. **Dimension** — the broad kind of information being measured.
2. **Category** — a coherent group of related concepts within that dimension.
3. **Feature** — one versioned question assessed against the title and abstract.

The hierarchy lives in `configs/features.yaml`, so it can be inspected and edited alongside each feature's ID, definition, applicable organism groups, question version, evidence instructions, criteria, and review status. `configs/feature_profiles.yaml` remains the explicit ordered list of features to score for each organism group. The Python loader validates and flattens the hierarchy for scoring; the Explorer consumes the same hierarchy to group feature summaries.

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

`schema_version: 2` identifies the nested YAML structure. The catalogue's feature/question versions continue to identify scoring meaning. A structure-only move does not invalidate successful decisions; a changed feature meaning does. Domain profiles continue to control the applicable feature set and its order.
