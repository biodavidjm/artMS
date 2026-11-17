# Converts the Protein ID column of the evidence file selected by the user to mod-site-specific notation: `ProteinID` to `ProteinID_AAnumber` notation

It enables the modified-peptide specific quantification by converting
the Protein column of the evidence file selected by the user to an
`ProteinID_AAnumber`notation. In this way, each of the modified peptides
can be quantified independently across conditions.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WARNING: we have detected a version of MaxQuant (\>1.6.3.0) outputs a\`
"Modified sequence" column of the evidence file that has two important
changes for the annotation of phosphorylation:

- Uses `p` instead of `(ph)`

- The modified residue (i.e. `STY`) is the residue on the right of the
  `p`, instead of the residue to the left of `(ph)`, as usual. We have
  introduced a modification to detect and address this issue, but we
  advice the user to double check both the new evidence file with the
  introduce new notation and the `-mapping.txt` file and check that
  there are no NA values for the notation of phophopeptides.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## Usage

``` r
artmsProtein2SiteConversion(
  evidence_file,
  ref_proteome_file,
  column_name = c("Leading razor protein", "Leading proteins", "Proteins"),
  output_file,
  mod_type,
  overwrite_evidence = FALSE,
  verbose = TRUE
)
```

## Arguments

- evidence_file:

  (char) The evidence file name and location

- ref_proteome_file:

  (char) The reference proteome used as database to search the
  `evidence.txt` file with MaxQuant. It will be used to map the modified
  peptide to the protein sequence and find the site location. Therefore,
  it does not use the MaxQuant's `Phospho (STY)Sites.txt`

- column_name:

  (char) The Protein Column Name to map. Options:

  - `Leadind razor protein` (default)

  - `Leading protein`

  - `Proteins` It only supports Uniprot Entry IDs and RefSeq, but it
    might work for other database IDs

- output_file:

  (char) Output file name (`ptmsites-evidence.txt` recommended)

- mod_type:

  (char) The posttranslational modification. Options:

  - `UB`: Protein Ubiquitination

  - `PH`: Protein Phosphorylation

  - `AC`: Protein Acetylation

  - `PTM:XXX:yy` : User defined PTM. Replace XXX with 1 or more 1-letter
    amino acid codes on which to find modifications (all uppercase).
    Replace yy with modification name used within the evidence file
    (require lowercase characters). Example: `PTM:STY:ph` will find
    modifications on aa S,T,Y with this format `_AAGGAPS(ph)PPPPVR_`.
    This would be equivalent to `mod_type = PH`

- overwrite_evidence:

  (logical) if \<output_file\> is the same as \<evidence_file\>,
  `overwrite_evidence = FALSE` (default) doesn't allow to overwrite the
  evidence file. Otherwise, `overwrite_evidence = TRUE` allows to
  overwrite the evidence_file (this option might be activated if the
  user allows to use the same `ptm-sites-evidence.txt` file to
  re-annotate all the Protein IDs columns)

- verbose:

  (logical) `TRUE` (default) shows function messages

## Value

(file) Return a new evidence file with the specified Protein id column
modified by adding the sequence site location(s) + postranslational
modification(s) to the uniprot entry / refseq id.

Output ID examples: `A34890_ph3`; `Q64890_ph24_ph456`;
`Q64890_ub34_ub129_ub234`; `Q64890_ac35`.

## Examples

``` r
# Testing warning if files are not submitted. 
artmsProtein2SiteConversion(evidence_file = NULL, ref_proteome_file = NULL, 
output_file = NULL)
#> [1] "Files must not be NULL"
```
