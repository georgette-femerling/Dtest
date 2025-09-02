import setuptools

with open("README.md", "r", encoding="utf-8") as fhand:
    long_description = fhand.read()

setuptools.setup(
    name="Dtest",
    version="0.1.0",
    author="Georgette Femerling",
    author_email="maria.femerling@mail.mcgill.com",
    description=("A command-line tool wrapper for moments.LD "
                "to test demographic models."),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/georgette-femerling/",
    project_urls={
        "Bug Tracker": "https://github.com/georgette-femerling/Dtest/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(include=['Dtest']),
    python_requires=">=3.6",
    entry_points={"console_scripts": ["dtest = Dtest.cli:main",]}
)