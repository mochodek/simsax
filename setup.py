from setuptools import setup

setup(name='simsax',
      version='0.1.1',
      description='Python SimSAX',
      url='https://github.com/mochodek/simsax',
      author='',
      author_email='',
      license='Apache-2.0',
      packages=[
        'simsax', 
        'simsax.alignment', 
        'sax'],
      install_requires=[
          'pandas',
          'numpy',
          'joblib',
          'biopython',
          'matplotlib'
      ],
      scripts=[
          'bin/simsax'],
      zip_safe=False)
