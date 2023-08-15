#!/bin/bash

R CMD BATCH 1a-namibia-incidence_preprocess.R &
R CMD BATCH 1b-namibia-incidence_preprocess-sens-obs.R &
R CMD BATCH 1c-namibia-incidence_preprocess-sens-nooverlap-target.R &
R CMD BATCH 1d-namibia-incidence_preprocess-sens-nooverlap-spill.R &
R CMD BATCH 1e-namibia-incidence_preprocess-sens-spzone-2km.R &
R CMD BATCH 1f-namibia-incidence_preprocess-sens-spzone-3km.R &
