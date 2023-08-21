import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="frag_hop",
    version="0.0.1",
    author="Laura Malo",
    author_email="laumaro95@gmail.com",
    description="Fragment Hopping",
    url="https://github.com/laumalo/frag_hop",
    packages=[
        'frag_hop',
        'frag_hop/data',
        'frag_hop/data/parameters',
        'frag_hop/replacement',
        'frag_hop/utils'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
