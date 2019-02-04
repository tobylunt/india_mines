# Mine detection from satellite imagery in India

This directory contains my effort to transfer learn a resnet
originally trained to classify satellite imagery using the Planet
Kaggle dataset from the challenge ["understanding the Amazon from
space](https://www.kaggle.com/c/planet-understanding-the-amazon-from-space). I
am adapting this CNN to use an expanded training dataset from known
mines in India, using publicly available Sentinel imagery from the
European Space Agency.

## Why should we care?

Mining activity in faraway places has massive ramifications on our
financial, emotional, and moral wellbeing. Being able to identify
mines accurately in remote areas, using pictures from space, would be
a powerful asset for tracking environmental degradation, key assets
for violent non-state actors, and the welfare of the rural poor.

This is a tough problem, mostly due to a lack of data at adequate
resolution and the massive computational resources
required. Fortunately modern techniques are reducing the barriers to
this kind of analysis.