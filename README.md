# hello-Sketch
hello-Sketch: An attempt of image retrieval based on sketch</br>

## introduction
This topic, 'sketch based image retrieval', is my subject that I work on in the scientific research training that I participated in during my sophomore year(presumably from Janurary 2019 to December 2019). In this work, I try to design a  hand-craft local feature that can be effective for image classification and retrieval task.</br>

My work can be summarized by following steps (for one RGB image, it may experiecnce following operations):</br>
1. get grayed
2. get edge extracted</br>
3. edge fitted with straight line</br>
4. a series local feature (say with shape [m,1]) extracted with a kind of feature design to compose the image's global feature (say with shape [m,n])</br>
5. the image's glocbal feature is feed into a [fisher vector](https://www.vlfeat.org/overview/encodings.html) to be encoded to a one-dimension vector (say with shape [a, 1])</br>

***every RGB image in dataset should go through all 1~5 process***</br>
***every sketch image in dataset should go through 2~5 process***</br>

***while all image in dataset have got their feature vector, companied with their labels, following tasks can be done***</br>
1. for image classification task, feed the feature vactors and labels to train a svm model.
2. for image retrieval task, the similarity of two images can be caculated using euclidean distance or other distance norm, according to the distance, retrival task can be fullfilled.
