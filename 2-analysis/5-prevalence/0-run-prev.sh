#!/bin/bash

R CMD BATCH 6a-namibia-prev-analysis-tmle3.R &
R CMD BATCH 6b-namibia-prev-analysis-tmle3-alt.R &
R CMD BATCH 6c-namibia-prev-analysis-tmle3-hotspot.R &
R CMD BATCH 6d-namibia-prev-analysis-serology.R &
R CMD BATCH 6e-namibia-prev-analysis-tmle3-dist.R &