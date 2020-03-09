# Multimedia-Systems

### Description
This collection of scripts was created for the scope of an assignment in a Fuzzy Systems cource of the Aristotle University of Thessaloniki during the 2019-20 academic year. The purpose of these scripts is to demonstrate the implementation of the baseline sequential DCT based (Lossy) JPEG encoder and decoder as described by ISO-IEC-10918-1-1993 standard.

![JPEG Encoded BitStream](https://github.com/kosletr/Multimedia-Systems/blob/master/Pics/Screenshot_8.jpg)

### Implementation

Implementation of the JPEG standard, as well as examples of its functionality are provided. The following block diagrams describe visually the given code.

 Encoder Block Diagram     | Decoder Block Diagram
:-------------------------:|:-------------------------:
![](https://github.com/kosletr/Multimedia-Systems/blob/master/Pics/Screenshot_1.jpg) |  ![](https://github.com/kosletr/Multimedia-Systems/blob/master/Pics/Screenshot_2.jpg)

Examples of the encoder/decoder process.

   Example 1               |        Example 2 
:-------------------------:|:-------------------------:
![](https://github.com/kosletr/Multimedia-Systems/blob/master/results/results/results_2.jpg) |  ![](https://github.com/kosletr/Multimedia-Systems/blob/master/results/results/results_19.jpg)

There is a trade off between high Compression and high Image Quality.

   Example 1               |        Example 2 
:-------------------------:|:-------------------------:
![](https://github.com/kosletr/Multimedia-Systems/blob/master/results/resultsB/resultsB_10.jpg) |  ![](https://github.com/kosletr/Multimedia-Systems/blob/master/results/resultsB/resultsB_25.jpg)

Moreover, various statistical results and graphs, about compression ratios, image entropy, bitstream sizes, quantization parameters etc. is also provided.

   Results Example 1       |  Results Example 2 
:-------------------------:|:-------------------------:
![](https://github.com/kosletr/Multimedia-Systems/blob/master/results/resultsB/resultsB_14.jpg) |  ![](https://github.com/kosletr/Multimedia-Systems/blob/master/results/resultsB/resultsB_15.jpg)


### Setup
The provided code was created using MATLAB R2019a, however older software versions should work fine. All of the .m scripts provided, are  commented for higher readability and maintenance.
