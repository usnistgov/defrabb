export const meta = {
  name: 'param-optimization-design',
  description: 'Design parameter optimization framework with adversarial review',
  phases: [
    { title: 'Design', detail: 'Architect parameter optimization framework' },
    { title: 'Critique', detail: 'Adversarial review for over/under-engineering' },
    { title: 'Refine', detail: 'Incorporate critique and finalize design' },
  ],
}

const DESIGN_SCHEMA = {
  type: 'object',
  required: ['architecture', 'components', 'rationale'],
  properties: {
    architecture: { type: 'string' },
    components: {
      type: 'array',
      items: {
        type: 'object',
        required: ['name', 'purpose', 'files'],
        properties: {
          name: { type: 'string' },
          purpose: { type: 'string' },
          files: { type: 'array', items: { type: 'string' } },
          implementation_notes: { type: 'string' },
        }
      }
    },
    rationale: { type: 'string' },
    complexity_estimate: { type: 'string' },
  }
}

const CRITIQUE_SCHEMA = {
  type: 'object',
  required: ['over_engineered', 'under_engineered', 'missing', 'verdict'],
  properties: {
    over_engineered: { type: 'array', items: { type: 'string' } },
    under_engineered: { type: 'array', items: { type: 'string' } },
    missing: { type: 'array', items: { type: 'string' } },
    risks: { type: 'array', items: { type: 'string' } },
    verdict: { type: 'string', enum: ['approve', 'revise', 'reject'] },
    recommendations: { type: 'string' },
  }
}

phase('Design')
log('Designing parameter optimization framework...')

const design_prompt = `Design the parameter optimization framework for DeFrABB v0.023.

Context from research:
- Dipcall params are broken (defined but never accessed)
- VCF processing and truvari bench params work well as patterns
- Exclusion params need profiles
- Output reuse is CRITICAL for 5-genome sweeps (dipcall takes 3-6 hours)

Requirements:
1. Fix dipcall parameter profile lookup
2. Add parameter sweep generation (analyses table rows from declarative config)
3. Enhance compare_evaluations.py for v5q baseline + scoring
4. Design output reuse for variant calls (reuse expensive dipcall/PAV runs)
5. Support multi-genome workflow: HG002 optimization then apply to 4 other genomes

Design principles:
- Build on existing infrastructure (VCF profiles, truvari params work well)
- Minimal changes to core pipeline (backward compatible)
- Incremental value (each component useful independently)
- Clear separation: param definition to sweep generation to evaluation to comparison

Output architecture overview, component list with files, rationale, complexity estimate.
Be concrete with file paths. Avoid over-engineering.`

const design = await agent(design_prompt, { phase: 'Design', schema: DESIGN_SCHEMA })

phase('Critique')
log('Running adversarial critiques in parallel...')

const design_json = JSON.stringify(design, null, 2)

const over_eng_prompt = 'Adversarially critique this parameter optimization design for OVER-engineering.\n\nDesign:\n' + design_json + '\n\nYour role: Find where the design is too complex, adds unnecessary abstraction, or tries to solve future problems we do not have yet.\n\nLook for:\n- Premature abstraction (helper classes for 1 use case)\n- Feature creep (nice-to-haves disguised as requirements)\n- Over-generalization (supporting scenarios we will not encounter)\n- Unnecessary layers (data passes through too many transformations)\n\nBe harsh but fair. Output specific components or decisions that are over-engineered.'

const under_eng_prompt = 'Adversarially critique this parameter optimization design for UNDER-engineering.\n\nDesign:\n' + design_json + '\n\nYour role: Find where the design is too simple, will break under real use, or ignores important edge cases.\n\nLook for:\n- Missing validation (schema, type checking, bounds)\n- No error handling (what breaks when inputs are wrong)\n- Scalability issues (works for 5 genomes, breaks at 50)\n- Missing observability (cannot debug when it fails)\n- Incomplete features (80 percent solution that forces manual workarounds)\n\nBe harsh but fair. Output specific gaps that will cause problems.'

const critiques = await parallel([
  () => agent(over_eng_prompt, { phase: 'Critique', schema: CRITIQUE_SCHEMA }),
  () => agent(under_eng_prompt, { phase: 'Critique', schema: CRITIQUE_SCHEMA }),
])

phase('Refine')
log('Refining design based on critiques...')

const refine_prompt = 'Refine the parameter optimization design based on adversarial critiques.\n\nOriginal design:\n' + design_json + '\n\nOver-engineering critique:\n' + JSON.stringify(critiques[0], null, 2) + '\n\nUnder-engineering critique:\n' + JSON.stringify(critiques[1], null, 2) + '\n\nProduce final design that:\n- Removes over-engineered components (simplify)\n- Adds under-engineered safeguards (validate and error-handle)\n- Balances complexity vs robustness\n- Stays focused on v0.023 goals (5-genome param optimization)\n\nOutput revised architecture and components with critique-driven changes highlighted.'

const final_design = await agent(refine_prompt, { phase: 'Refine', schema: DESIGN_SCHEMA })

return {
  initial_design: design,
  over_engineering_critique: critiques[0],
  under_engineering_critique: critiques[1],
  final_design: final_design,
}
