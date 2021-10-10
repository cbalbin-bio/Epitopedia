from setuptools import setup, find_packages


setup(
    name="Epitopedia",
    version="1.1",
    description="Structural molecular mimicry discovery pipeline",
    author="Christian Balbin",
    author_email="cbalbin@fiu.edu",
    packages=find_packages(),
    python_requires=">=3.9, <4",
    install_requires=["flask", "gemmi", "biopython", "rich", "dataclasses-json"],
    entry_points={  # Optional
        "console_scripts": [
            "generate_database=epitopedia.generate_database:main",
            "run_epitopedia=epitopedia.run_epitopedia:main",
        ],
    },
)
