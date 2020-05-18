[![License](https://img.shields.io/badge/License-Apache%202.0-yellowgreen.svg)](https://opensource.org/licenses/Apache-2.0)  

# A Probabilistic Calibration Framework for the Abaqus Superelastic Material Model


## Introduction

![Probabilistic Framework Abstract](/assets/graphical_abstract.png)

We implement a Bayesian Inference approach to calibrate the material parameters of a model for the superelastic deformation of NiTi shape memory alloy. We specify a diamond-shaped specimen geometry that is suited to calibrate both tensile and compressive material parameters from a single test. We adopt the Bayesian Inference calibration scheme to take full-field strain measurements obtained using digital image correlation together with global load data as an input for calibration. The calibration itself is performed by comparing the experimentally measured quantities of interest -- strain data at select locations and global load -- with the corresponding results from a simulation library. We present a machine learning based approach to enrich the simulation library and improve the calibration accuracy. This approach is versatile and can be used to calibrate other models of superelastic deformation from data obtained using various modalities. This probabilistic calibration approach can become an integral part of a framework to assess and communicate the risk-informed credibility of simulations performed in the design of superelastic NiTi articles such as medical devices.

Matlab code to execute the full calibration pipeline is presented here. For details see the research article describing the method.

### Publication

Harshad M. Paranjape, Kenneth I. Aycock, Craig Bonsignore, Jason D. Weaver, Brent A. Craven, Thomas W. Duerig. "A Probabilistic Approach with Built-in Uncertainty Quantification for the Calibration of a Superelastic Constitutive Model from Full-field Strain Data." Under review.

### Prerequisites

Below are pre-requisite tools to perform the material parameter calibration described here.
* Matlab with Image Processing Toolbox and Statistics and Machine Learning Toolbox.
* A diamond specimen with the geometry specified in the engineering drawing `data/NDC-56-03018_rev3_25102018_00.pdf`. This specimen should be laser-cut and heat treated from the material to which the model is to be calibrated.
* Instron or similar load frame to perform tension experiments on the specimen.
* A 2D digital image correlation (DIC) setup including tools for speckle pattern application on the diamond, digital camera, lighting, and other supporting equipment.
* NCORR DIC data processing software for Matlab.
* Matlab source code in this repository.



This is an example of how to list things you need to use the software and how to install them.
* npm
```sh
npm install npm@latest -g
```

### Installation

1. Get a free API Key at [https://example.com](https://example.com)
2. Clone the repo
```sh
git clone https://github.com/your_username_/Project-Name.git
```
3. Install NPM packages
```sh
npm install
```
4. Enter your API in `config.js`
```JS
const API_KEY = 'ENTER YOUR API';
```



<!-- USAGE EXAMPLES -->
## Usage

Use this space to show useful examples of how a project can be used. Additional screenshots, code examples and demos work well in this space. You may also link to more resources.

_For more examples, please refer to the [Documentation](https://example.com)_



<!-- ROADMAP -->
## Roadmap

See the [open issues](https://github.com/othneildrew/Best-README-Template/issues) for a list of proposed features (and known issues).



<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to be learn, inspire, and create. Any contributions you make are **greatly appreciated**.

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.



<!-- CONTACT -->
## Contact

Your Name - [@your_twitter](https://twitter.com/your_username) - email@example.com

Project Link: [https://github.com/your_username/repo_name](https://github.com/your_username/repo_name)



<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements
* [GitHub Emoji Cheat Sheet](https://www.webpagefx.com/tools/emoji-cheat-sheet)
* [Img Shields](https://shields.io)
* [Choose an Open Source License](https://choosealicense.com)
* [GitHub Pages](https://pages.github.com)
* [Animate.css](https://daneden.github.io/animate.css)
* [Loaders.css](https://connoratherton.com/loaders)
* [Slick Carousel](https://kenwheeler.github.io/slick)
* [Smooth Scroll](https://github.com/cferdinandi/smooth-scroll)
* [Sticky Kit](http://leafo.net/sticky-kit)
* [JVectorMap](http://jvectormap.com)
* [Font Awesome](https://fontawesome.com)





<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/othneildrew/Best-README-Template.svg?style=flat-square
[contributors-url]: https://github.com/othneildrew/Best-README-Template/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/othneildrew/Best-README-Template.svg?style=flat-square
[forks-url]: https://github.com/othneildrew/Best-README-Template/network/members
[stars-shield]: https://img.shields.io/github/stars/othneildrew/Best-README-Template.svg?style=flat-square
[stars-url]: https://github.com/othneildrew/Best-README-Template/stargazers
[issues-shield]: https://img.shields.io/github/issues/othneildrew/Best-README-Template.svg?style=flat-square
[issues-url]: https://github.com/othneildrew/Best-README-Template/issues
[license-shield]: https://img.shields.io/github/license/othneildrew/Best-README-Template.svg?style=flat-square
[license-url]: https://github.com/othneildrew/Best-README-Template/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=flat-square&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/othneildrew
[product-screenshot]: images/screenshot.png
