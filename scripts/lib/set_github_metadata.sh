#!/usr/bin/env bash
# Set the GitHub description and topics once the repository has a remote.
# Not run by the pipeline; run it after `git push`.
#
# Usage: scripts/lib/set_github_metadata.sh [owner/repo]
set -euo pipefail
REPO="${1:-}"
ARGS=()
if [[ -n "$REPO" ]]; then ARGS=("$REPO"); fi
DESCRIPTION="A pre-registered reanalysis of the MetaboLights MTBLS9264 stool metabolomics of people living with HIV and healthy controls, processed from the raw Orbitrap spectra with MassCube one file at a time within a 10 GB disk budget, with an antiretroviral positive control, QC-based filtering, limma tests and a claim-by-claim replication of the published findings."
TOPICS=(metabolomics untargeted-metabolomics lc-ms masscube hiv gut-microbiome stool-metabolomics metabolights
        reproducibility pre-registration replication limma ropls opls-da kegg-pathways quality-control
        thermorawfileparser orbitrap bioinformatics r python)
gh repo edit "${ARGS[@]+"${ARGS[@]}"}" --description "$DESCRIPTION"
for t in "${TOPICS[@]}"; do
    gh repo edit "${ARGS[@]+"${ARGS[@]}"}" --add-topic "$t"
done
echo "description and ${#TOPICS[@]} topics set"
