
Only the **positions** section of the JSON is used.  
Each position entry corresponds to a variant present in the unannotated VCF.

---

## Transcript Selection

### Primary method
- Use the latest **MANE Select** RefSeq transcripts  
- Validate Nirvana’s `isManeSelect` flag against the MANE GTF

### Fallback logic
Used when MANE is not available:

- Prefer transcripts where:  
  - `biotype == "mRNA"`  
  - `isCanonical == true`  

This ensures a single, deterministic transcript is selected for each variant.

---

## Output

A TSV file where each row represents a variant.

### Variant fields
- Chromosome  
- Position  
- refAllele  
- altAllele(s)  
- Filters  
- mappingQuality  
- cytogeneticBand  
- genotype  
- variantFrequencies  
- totalDepth  
- alleleDepths  
- somaticQuality  
- variantType  
- hgvsg  

### ClinVar fields
- ID  
- significance  
- reviewStatus  
- Only entries matching the variant ref/alt are included  
- Multiple entries are combined with a delimiter  

### Population frequencies
- dbsnp ID  
- globalMinorAlleleFrequency (if matching altAllele)  
- gnomad allAF  
- oneKg allAF  
- topmed allAF  

### COSMIC
- COSMIC IDs (if any)

### Transcript fields
- Transcript name (RefSeq)  
- biotype  
- exon or intron notation  
- hgnc  
- consequence  
- hgvsc  
- hgvsp (if present)

---

## Validation

Validation should include:

- Unit tests for parsing and transcript selection  
- Simple and complex representative variants  
- Manual spot-checking to confirm correctness
