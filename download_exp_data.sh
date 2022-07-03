#/bin/bash

for i in $(seq 1 3)
do
    PATH_TEST="~/experiments/test$i/GRIEF-DE"
    echo $PATH_TEST 
    cd experiments
    mkdir test$i
    cd ..
    # scp -r -i ~/.ssh/gcloud cicero_samuel_mendes@35.238.165.151:"$PATH_TEST/results $PATH_TEST/test_description.txt $PATH_TEST/run.sh" ~/Dev/msc/GRIEF-DE-c++/experiments/test$i
    scp -r -i ~/.ssh/gcloud cicero_samuel_mendes@35.238.165.151:"$PATH_TEST/results" ~/Dev/msc/GRIEF-DE-c++/experiments/test$i
done