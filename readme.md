Semantic Soft Segmentation
=================================================================

This repository includes the spectral segmentation approach presented in

    Yagiz Aksoy, Tae-Hyun Oh, Sylvain Paris, Marc Pollefeys and Wojciech Matusik, "Semantic Soft Segmentation", ACM Transactions on Graphics (Proc. SIGGRAPH), 2018 

The network for semantic feature generation can be found [[here](https://github.com/iyah4888/SIGGRAPH18SSS)].

Please refer to the [[project page](http://people.inf.ethz.ch/aksoyy/sss/)] for more information and example data.

The Superpixels class requires [[ImageGraphs](http://www.mathworks.com/matlabcentral/fileexchange/53614-image-graphs)].

License and citation
------------

This toolbox is provided for academic use only.
If you use this code, please cite our paper:

    @ARTICLE{sss,
    author={Ya\u{g}{\i}z Aksoy and Tae-Hyun Oh and Sylvain Paris and Marc Pollefeys and Wojciech Matusik},
    title={Semantic Soft Segmentation},
    journal={ACM Transactions on Graphics (Proc. SIGGRAPH)},
    year={2018},
    pages = {72:1-72:13},
    volume = {37},
    number = {4}
    }

Credit
------------
Parts of this implementation are taken from

    Anat Levin, Alex Rav-Acha, Dani Lischinski, "Spectral Matting", IEEE TPAMI, 2008

The original source code for Spectral Matting can be found [[here](http://www.vision.huji.ac.il/SpectralMatting/)].