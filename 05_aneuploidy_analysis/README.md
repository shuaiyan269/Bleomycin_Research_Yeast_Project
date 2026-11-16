# 05_aneuploidy_analysis
## Server Environment (Bash)

Step 1:
For the post-mapping depth data `all_samples.coverage_summary.txt`, run `summary_median_depth_normalized.sh` to calculate the copy number, defined as (chromosome median depth / median depth) × 2.

Step 2:
For the post-mapping per-sample depth files `sample.depth.txt.gz` (not provided here due to their large size), run `depth_sld_win.pl` to obtain the 10 kb sliding-window depth data for each sample.
