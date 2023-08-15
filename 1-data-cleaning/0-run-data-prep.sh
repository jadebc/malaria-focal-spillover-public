#!/bin/bash

R CMD BATCH 5-clean-namibia-xs.R &

R CMD BATCH 6a-data-prep.R &
R CMD BATCH 6b-data-prep-short.R &
R CMD BATCH 6c-data-prep-long.R &

R CMD BATCH 7a-data-prep-short-sens-obs.R &
R CMD BATCH 7b-data-prep-long-sens-obs.R &
R CMD BATCH 7c-data-prep-short-sens-spzone-2km.R &
R CMD BATCH 7d-data-prep-long-sens-spzone-2km.R &
R CMD BATCH 7e-data-prep-short-sens-spzone-3km.R &
R CMD BATCH 7f-data-prep-long-sens-spzone-3km.R &
R CMD BATCH 8-data-prep-sens-no-overlap.R &
