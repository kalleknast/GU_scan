# GU_scan
Script and functions for in vitro mapping of effective neural connectivity using UV controlled glutamate uncaging (GU).
Efficiently samples neural connections, combining the opposing requirements for high spatial resolution, large coverage and quick execution time.

## Background to GU
The scripts were developed for an in vitro neurophysiology experiment done in order to map connections between neurons within the bed nucleus of the stria terminalis (BNST) using patch recordings and ultraviolet (UV) glutamate uncaging (GU). This works through illuminating a circular region of the BNST with UV-light while recording the membrane potential in one cell. UV-light releases glutamate which in turn activates the neurons in the illuminated region. If an activated cell is connected to the recorded cell then this is detected as changes in the recorded cell's membrane potential. Through sequentially illuminating several regions one can map the connections between the a given cell and those regions. 

### Reference
Turesson HK, Rodríguez-Sierra OE, Pare D. Intrinsic connections in the anterior part of the bed nucleus of the stria terminalis. Journal of Neurophysiology. 2013;109(10):2438-2450.<br/>
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3653047/

## The problem
The spatial resolution of this mapping is controlled by the size of the region illuminated. Smaller region, higher resolution. The goal is often to map as big area as possible with as high resolution as possible. However, the time that a neuron can be recorded is limited and each region sampled takes time. This limits the size and resolution of the sampled area.

## Solution
In order to maximize the sampled region and resolution, I developed a method similar many compression methods. Starting with max size of illuminated regions the area of interest in sampled and connected regions registered. With a 50% smaller illumination region the regions that in the previous step were registered as connected are sampled. This is repeated until the smallest region size, and max resolution is reached.
