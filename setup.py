from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='emp_hagai',
    version="0.1",
    author="Hagai Levi",
    author_email="hagai.levi.007@gmail.com",
    long_description=long_description,
    long_description_content_type="text/markdown",
    description='EMP: EMpirical Pipeline for correcting enrichment scores of network-based module discovery algorithms',
    url='https://github.com/hag007/DOMINO',
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Linux",
    ],
    packages = find_packages(),
    install_requires=[
        'fastsemsim==1.0.0',
        'goatools==0.9.9',
        'numpy==1.18.1',
        'pandas==1.0.1',
        'pcst-fast==1.0.7',
        'scikit-learn==0.22.1',
        'scipy==1.4.1',
        'simplejson==3.17.0',
        'statsmodels==0.11.0',
        'wget==3.2'
    ],
    entry_points = {
        "console_scripts": [
            "generate_permuted_solutions=src.emp.generate_permuted_solutions:main",
            "aggregate_bg_distribution=src.emp.aggregate_bg_distribution:main",
            "add_go_metadata=src.emp.add_go_metadata:main",
            "calculate_significance=src.emp.calculate_significance:main",
            "bg_dist_to_pval=src.emp.bg_dist_to_pval:main",
        ]
    }



)