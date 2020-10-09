import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pykasso",
    version="0.1.0",
    author="Francois Miville",
    author_email="francois.miville@unine.ch",
	url="https://github.com/randlab/pyKasso",
    description="Python project intended to simulate stochastic karst network.",
    long_description=long_description,
	long_description_content_type="text/markdown",
	classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    keywords="simulation karstic stochastic",
    packages=setuptools.find_packages(),
    python_requires='>=3.6',
    include_package_data=True
)
