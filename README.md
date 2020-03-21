# Julia-Package-convalescent_sera.jl
# This notebook summarizes initial thoughts about the "convalescent sera" project

Initial thoughts by Shane, Marijn, Sarat, Yoni, Peter  + (required: public health experts, blood bank experts, actual clinical data)

Last updated March 19, 2020.

**Background**: There is an option of blood transfusions (convalescent sera) from recovered patients to sick or highly sucscuipulte ones. We want to quickly quantify the effect of this and help the medical community via mathematical models that will help drive decisions and predict the future for policy makers.

See [this paper](https://annals.org/aim/fullarticle/729754/meta-analysis-convalescent-blood-products-spanish-influenza-pneumonia-future-h5n1) and this [public release](https://hub.jhu.edu/2020/03/13/covid-19-antibody-sera-arturo-casadevall/) summarizing [this paper](https://www.jci.org/articles/view/138003).

### Who will use our type of results?

1. Policy makers (good because they can make better informed decisions).
1. Clinicians (tradeoff between this type of treatment and others... at the patient and public health levels).
1. Knoweledge for general public (knowing where things are going and why you need to donate your blood is important).
1. Blood bank managers (?)

### Here are some questions that we may want to try to answer:

1. At the patient (or family, or hospital, or group) level: how is it best to apply convalesecent sera, given the positive/negative tradeoffs, costs, benefit, and risks that will be experienced in testing.

1. To which groups should this be applied first (or in what priority). More generally, what are sensible application protocols.

1. Assuming the convalescent sera is a somewhat viable option, how will it help to "flatten the curve". (quantitativly).

It appears that without further knoweldge of parameters, setting out how to answer the third question is the first thing to consider.

### Projects (or sub-projects)

1. Patient level model.
1. Population based decision trees.
1. Single compartment level (SEIR model): A good basic read if you don't know SEIR is [this blog post](https://triplebyte.com/blog/modeling-infectious-diseases)
1. National level models
1. Global level models

# What questions/parameters would we need from the medical field?

(This is based on knowledge that will probably be aquired during clinical testing in weeks/months to come)

1. Are the "recovered" all in the same pool (except for blood type) or do the anti-bodies "fade-out"... and if so how...
1. What are properties of CS? What are the risks? What are the attributes of success or failure? 
1. Is it sensible for individuals to recieve multiple doses - especially for prevention if dealing with infected often (e.g. medical staff)? 

