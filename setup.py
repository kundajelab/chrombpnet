from setuptools import setup,find_packages

#generate install_requires from requirements.txt file
install_requires=open('requirements.txt','r').read().strip().split('\n')
print(f"install_requires:{install_requires}")


config = {
    'name': 'chrombpnet',
    'author': 'Kundaje lab',
    'author_email': 'anusri @ stanford.edu',
    'license': 'MIT',
    'include_package_data': True,
    'description': 'chrombpnet predicts chromatin accessibility from sequence',
    'download_url': 'https://github.com/kundajelab/chrombpnet',
    'version': '0.1',
    'packages': find_packages(),
    'python_requires': '>=3.6',
    'install_requires': install_requires,
    'zip_safe': False,
    'scripts':['chrombpnet/helpers/make_gc_matched_negatives/make_gc_matched_negatives.sh',
               'chrombpnet/helpers/preprocessing/bam_to_bigwig.sh',
               'chrombpnet/training/models/bpnet_model.py',
               'chrombpnet/training/models/chrombpnet_with_bias_model.py',
               'step1_download_bams_and_peaks.sh',
               'step2_make_bigwigs_from_bams.sh',
               'step3_get_background_regions.sh',
               'step4_train_bias_model.sh',
               'step5_interpret_bias_model.sh',
               'step6_train_chrombpnet_model.sh',
               'step7_interpret_chrombpnet_model.sh'
    ],
    'entry_points': {'console_scripts': [
        'chrombpnet_genomewide_gc = chrombpnet.helpers.make_gc_matched_negatives.get_genomewide_gc_buckets.get_genomewide_gc_bins:main',
        'chrombpnet_gc_content_foreground = chrombpnet.helpers.make_gc_matched_negatives.get_gc_content:main',
        'chrombpnet_gc_matched_negatives = chrombpnet.helpers.make_gc_matched_negatives.get_gc_matched_negatives:main',
        'chrombpnet_pwm_from_bigwig = chrombpnet.helpers.preprocessing.analysis.build_pwm_from_bigwig:main',        
        'chrombpnet_make_splits = chrombpnet.helpers.make_chr_splits.splits:main',
        'chrombpnet_hyperparams = chrombpnet.helpers.hyperparameters.find_chrombpnet_hyperparams:main',
        'chrombpnet_bias_hyperparams = chrombpnet.helpers.hyperparameters.find_bias_hyperparams:main',
        'chrombpnet_train = chrombpnet.training.train:main',
        'chrombpnet_predict = chrombpnet.training.predict:main',
        'chrombpnet_metrics = chrombpnet.training.metrics:main',
        'chrombpnet_deepshap = chrombpnet.evaluation.interpret.interpret:main',
        'chrombpnet_modisco = chrombpnet.evaluation.modisco.run_modisco:main',
        'chrombpnet_marginal_footprints = chrombpnet.evaluation.marginal_footprints.marginal_footprinting:main',
        'chrombpnet_tomtom_hits = chrombpnet.evaluation.modisco.fetch_tomtom:main',
        'chrombpnet_score_snps = chrombpnet.evaluation.variant_effect_prediction.snp_scoring:main',
        'chrombpnet_srcdir = chrombpnet.get_package_dir:main']}}

if __name__== '__main__':
    setup(**config)
