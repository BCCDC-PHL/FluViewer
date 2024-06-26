from setuptools import setup, find_namespace_packages

import fluviewer

setup(
    name='fluviewer',
    version=fluviewer.__version__,
    packages=find_namespace_packages(),
    entry_points={
        "console_scripts": [
            "FluViewer = fluviewer.fluviewer:main",
            "fluviewer = fluviewer.fluviewer:main",
        ]
    },
    scripts=[],
    package_data={
    },
    python_requires='>=3.7',
    install_requires=[
    ],
    description='',
    url='https://github.com/BCCDC-PHL/FluViewer',
    author='Kevin Kuchinski',
    author_email='',
    include_package_data=True,
    keywords=[],
    zip_safe=False
)
