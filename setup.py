import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="fasano-franceschini-test",
    version="1.0",
    author="Ian Chow",
    author_email = "ichow9@uwo.ca",
    description = "A Python implementation of the multivariate two-sample Kolmogorov-Smirnov test described by Fasano and Franceschini (1987)",
    long_description = long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    packages=["fftest"],
    setup_requires=['numpy>=1.8'],
    install_requires=['numpy']
    )