---
name: viroflash Evidence Ledger
description: An offline evidence-first review system for viral candidate reports.
colors:
  evidence-violet: "#6550a5"
  evidence-mint: "#1c775f"
  caution-amber: "#8a5b12"
  review-rose: "#a63d59"
  paper: "#f7f8f5"
  surface: "#ffffff"
  ink: "#182029"
  muted-ink: "#59646f"
  rule: "#d8ddd9"
  strong-rule: "#aeb8b3"
typography:
  display:
    fontFamily: "ui-sans-serif, -apple-system, BlinkMacSystemFont, Segoe UI, sans-serif"
    fontSize: "clamp(2.2rem, 4vw, 3.8rem)"
    fontWeight: 760
    lineHeight: 0.98
    letterSpacing: "-0.035em"
  headline:
    fontFamily: "ui-sans-serif, -apple-system, BlinkMacSystemFont, Segoe UI, sans-serif"
    fontSize: "1.45rem"
    fontWeight: 700
    lineHeight: 1.2
    letterSpacing: "-0.02em"
  body:
    fontFamily: "ui-sans-serif, -apple-system, BlinkMacSystemFont, Segoe UI, sans-serif"
    fontSize: "1rem"
    fontWeight: 400
    lineHeight: 1.55
  label:
    fontFamily: "ui-sans-serif, -apple-system, BlinkMacSystemFont, Segoe UI, sans-serif"
    fontSize: "0.75rem"
    fontWeight: 720
    lineHeight: 1.35
rounded:
  control: "10px"
  action: "12px"
  status: "999px"
spacing:
  compact: "0.5rem"
  control: "0.75rem"
  section: "2rem"
  major: "4.5rem"
components:
  button-primary:
    backgroundColor: "{colors.ink}"
    textColor: "{colors.surface}"
    rounded: "{rounded.action}"
    padding: "0.72rem 1rem"
  input:
    backgroundColor: "{colors.surface}"
    textColor: "{colors.ink}"
    rounded: "{rounded.control}"
    padding: "0.65rem 0.75rem"
  status-pass:
    backgroundColor: "#dff2e9"
    textColor: "#105c48"
    rounded: "{rounded.status}"
    padding: "0.28rem 0.6rem"
---

# Design System: viroflash Evidence Ledger

## Overview

**Creative North Star: "The Evidence Ledger"**

The system behaves like a precise review worksheet rather than a diagnostic dashboard. Paper-white fields, slate text, hairline rules, and narrow evidence bands keep dense scientific content calm while preserving unmistakable status and hierarchy. Interaction supports inspection; it never substitutes decoration for evidence.

The visual language is intentionally restrained because candidate uncertainty, QC limitations, and unresolved reference groups must remain legible beside strong numerical evidence.

**Key Characteristics:**

- Evidence and interpretation boundaries share the same visual level.
- Dense information is organized by rules, alignment, and rhythm rather than nested cards.
- Mint, rose, amber, and violet are semantic accents; every state also has a text label.
- The same report remains useful on desktop, mobile, and paper.

## Colors

The palette keeps the reading field achromatic and confines saturation to evidence, state, focus, and action.

### Primary

- **Evidence Violet:** Primary links, disclosure controls, and quantitative breadth rails.
- **Evidence Mint:** PASS state and distributed-evidence marks.

### Secondary

- **Caution Amber:** Unevaluated QC and below-threshold states.
- **Review Rose:** Not-significant states and statistical review signals.

### Neutral

- **Paper:** The page field used behind every report section.
- **Surface:** Candidate and control surfaces.
- **Ledger Ink:** Primary text and primary action fill.
- **Muted Ink:** Definitions, units, and supporting context.
- **Rule / Strong Rule:** Data-row separators and major section boundaries.

**The Narrow Spectrum Rule.** Saturated colors appear as narrow bands, short rails, status marks, or small chips; they never flood a scientific content region.

**The Redundant State Rule.** Color never carries decision meaning alone. Pair every semantic color with an explicit label and stable position.

## Typography

**Display Font:** Native UI sans stack
**Body Font:** Native UI sans stack
**Label/Mono Font:** Native UI sans for labels; SFMono/Consolas fallback stack for contract identifiers

**Character:** Compact, direct, and instrument-like. Tabular numerals stabilize scientific values while ordinary prose remains easy to scan.

### Hierarchy

- **Display** (760, responsive 2.2–3.8rem, 0.98): Sample identity only.
- **Headline** (700, 1.45rem, 1.2): Major review regions.
- **Title** (700, approximately 1rem–1.35rem): Candidates and audit groups.
- **Body** (400, 1rem, 1.55): Interpretation text with a 65–72 character measure.
- **Label** (720, 0.75rem): Controls, status, metrics, and field names.

**The Unit-Beside-Value Rule.** Units, denominators, and threshold semantics stay adjacent to the numerical value they constrain.

## Layout

The report uses a wide ledger at a maximum width of 1480px. Major regions are separated by hairline rules. Desktop layouts use two-column orientation, six-column run metadata, and four-column evidence/audit groups. At 1000px these collapse to two or three columns; at 640px they become a single reading flow, while short metric groups retain a two-column grid.

Spacing separates interpretation phases: compact spacing inside facts, approximately 2rem inside sections, and 3–4.5rem before major review regions. Responsive behavior changes structure rather than shrinking type continuously.

## Elevation & Depth

The system is flat by default. Borders, tonal fields, and section rhythm establish depth. Buttons may lift by one pixel with a diffuse ambient shadow on hover; scientific evidence containers do not float.

**The Flat Evidence Rule.** Data surfaces remain flat and aligned so visual elevation cannot be mistaken for evidence strength.

## Shapes

Controls use gently curved 10–12px corners. Small statuses use capsule geometry because they are compact labels, not containers. Evidence rails and section rules stay square and linear to preserve the ledger character.

## Components

### Buttons

- **Primary:** Ledger Ink surface, white text, 12px radius, compact horizontal padding.
- **Secondary:** Transparent surface with a strong neutral rule.
- **Hover / Focus:** One-pixel lift for pointer hover; a visible 3px blue focus outline for keyboard users.

### Chips

- **Style:** Compact uppercase text with either a semantic tonal fill or a neutral one-pixel outline.
- **State:** Always includes readable status text; no color-only variants.

### Cards / Containers

- **Corner Style:** Square scientific evidence regions; cards are not the page scaffold.
- **Background:** White surface on paper field.
- **Shadow Strategy:** None.
- **Border:** Horizontal rules and a narrow horizontal state band where candidate state needs rapid scanning.

### Inputs / Fields

- **Style:** White fill, strong neutral stroke, 10px radius, direct labels above controls.
- **Focus:** Shared 3px blue focus outline.

### Evidence Rails

- **Breadth:** Violet proportional fill with an explicit threshold marker and exact percentage label.
- **Distributed windows:** Ten square segments used only as a quantity indicator, paired with text that states positions are not encoded.

## Do's and Don'ts

### Do:

- **Do** put uncertainty and interpretation boundaries beside the evidence they constrain.
- **Do** use tabular numerals and explicit denominators for scientific values.
- **Do** preserve text labels, keyboard focus, responsive reflow, and print behavior.
- **Do** use the horizontal evidence band as a small scan aid, never as the primary decision signal.

### Don't:

- **Don't** present a candidate gate as a sample-level or clinical conclusion.
- **Don't** turn OR-group members into separate positive calls.
- **Don't** use large saturated regions, gradients in text, decorative glass, or floating evidence cards.
- **Don't** use vertical colored side stripes to mark candidates.
