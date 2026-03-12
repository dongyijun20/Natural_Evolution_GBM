Dear Reviewer,

Thank you very much for your insightful comment regarding potential pseudoreplication. We fully agree that treating individual cells as independent observations can inflate significance when data are hierarchically structured within patients.

In response, we have revised the statistical framework used in the Figure 1-5 analysis scripts and now retain a linear mixed-effects model (LMM) as the only inferential method for cross-group p-value calculation.

### What was changed

1. **Patient-aware modeling**
   - We model lesion/group effects using LMMs with patient as a random intercept:
   - `outcome ~ lesion + (1 | patient)`
   - This directly accounts for within-patient correlation and avoids pseudoreplication.

2. **Unit of replication**
   - For cross-group comparisons, inference is based on patient-level replication structure rather than treating each cell as an independent biological replicate.
   - The tested fixed effect is the lesion/group term, while patient-level heterogeneity is captured by the random effect.

3. **Single p-value method**
   - We removed parallel non-mixed tests (for example, paired Wilcoxon and paired t-tests) from the Figure 1-5 plotting workflow where these comparisons were previously reported together.
   - Figure captions now report only LMM-derived p-values.

### Why this addresses the concern

This update aligns the degrees of freedom with the number of biological subjects and yields more conservative, statistically rigorous inference. Importantly, while p-values are generally more conservative after this revision, the major biological trends and conclusions remain consistent with the original manuscript.

We sincerely appreciate this suggestion, which has strengthened the statistical rigor and interpretability of our analyses.
