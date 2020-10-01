# Thesis
육각그리드 IL-SOM을 이용한 개체시각화
---------------


자기조직화지도(SOM)는 고차원 다변량 자료를 저차원 그리드 공간에 축약하여 시각화하는 비지도학습 신경망 모형의 일종이다(Kohonen, 1998).SOM은 개체를 각 승자노드의 중심에 할당하는데 개체의 연속성을 보존한다는 관점에서는 그리드 공간에 승자노드를 중심으로 개체들이 적절히 퍼지도록 분포시키는 게 바람직하다. Um(2003)은 개체벡터의 승자노드와 그 인접노드의 중량벡터에 대한 가능도를 이용하여 그리드 공간상에 개체벡터의 상대적인 위치를 표현하도록 SOM을 확장한 IL-SOM(Interpolating Likelihood for SOM)을 제안하였다. SOM에서는 사각형 또는 육각형 그리드를 주로 사용하는데 UM(2003)은 이차원 그리드 공간에 사각그리드와 가상노드를 통해 확장한 확장형 사각그리드를 구현하였고 개체표현지수(Um, 2003)가 최소화되는 그리드의 크기를 선택하면 개체를 가장 잘 시각화한다는 결론을 내렸다. 본 연구에서는 이차원 그리드 공간에 육각그리드와 가상노드를 통해 확장한 확장형 육각그리드를 사용하는 IL-SOM을 구현하였다. 모의실험을 통해 사각그리드, 확장형 사각그리드, 육각그리드와 확장형 육각그리드의 개체표현지수 및 자료의 시각화 결과를 비교하였다. 또한 이 방법들을 4개의 그룹으로 구성된 게 자료에 적용하여 개체표현지수를 최소화하는 그리드의 종류를 확인하고 이를 통해 자료를 시각화한 결과를 비교하였다. 사각그리드와 육각그리드 모두 가상노드가 없는 경우보다 있는 경우가 개체표현지수가 작아 더 좋은 결과를 보였으며 육각그리드가 사각그리드보다 대체로 나은 결과를 보였다.

Visualization Using Hexagonal Grid IL-SOM
-----------------
  Self-organizing maps(SOM) is an unsupervised learning neural network model to visualize high dimensional data using a low dimensional grid space(Kohonen, 1990). SOM assigns objects to the center of the winning node, but it would be better to spread them all over the grid space to preserve continuity of the objects. Um(2003) proposed IL-SOM(Interpolating Likelihood for SOM) to express relative location of object vectors in a grid space using the likelihood of weight vectors of the winning node and its neighbor nodes. SOM usually uses the rectangular or the hexagonal grid. In a 2-dimensional rectangular grid, Um(2003) derived IL-SOM and extended it using a virtual nodes. In this study, IL-SOM is derived in a 2-dimensional hexagonal grid and extended hexagonal grid. These four methods are also compared with the crab data. In both rectangular and hexagonal grid, the extension using virtual nodes performed better. Performance of the hexagonal grid is usually better than the one of the rectangular grid.
  


## 참고 자료

* [논문 원문](https://drive.google.com/file/d/1MQvHF_roYARn-eJHJXq2UptAskAG-azw/view?usp=sharing)
* [KCI 등재 자료](https://www.dropbox.com/s/ld13n32o9diu500/%EC%9C%A1%EA%B0%81%EA%B7%B8%EB%A6%AC%EB%93%9C%20IL-SOM%EC%9D%84%20%EC%9D%B4%EC%9A%A9%ED%95%9C%20%EA%B0%9C%EC%B2%B4%20%EC%8B%9C%EA%B0%81%ED%99%94_KCI.pdf?dl=0)

