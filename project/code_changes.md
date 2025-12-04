## CODE CHANGES

lefse.py:
    - Removed R dependency and replaces R-based statistical tests with pure Python/SciPy implementations 
    - Adds a new test_cohens_d() function that replaces the R-based test_lda_r() function with bootstrapped    Cohen's d effect size calculation

lefse_plot_res.py:
    - Modified Cohen's D version with new features (effect size interpretation function, threshold lines, shaded regions showing effect size categories, 4-column output with effect labels) and different defaults (DPI 150, larger dimensions, title changes)

lefse_run.py:
    - Replaces LDA naming with Cohen's d
    - Adds better error handling (try-except blocks around feature testing, tracking of rejected features)
    - Removing verbose print statements for cleaner output
    - Cohen's d replaces LDA as the effect size measure
    - Update imports

AbundaceTable.py:
    - Update imports for new structure, import numpy to use numpy array, port from python2 to python3, use modern to_json method

ValidateData.py:
    - Fix syntax errors, port from python2 to python3
lefse_format_input.py, lefse_plot_cladogram.py:
    - Update imports

lefse_plot_features.py:
    - Fix UnboundLocalError

lefse2circlader.py & qiime2lefse.py:
    - No changes

bioconda-lefse_run.sh -> evaluate_mouse_data.sh
    - Modified to fit the mouse data
