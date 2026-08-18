# pixelAV_Generation

## Step 1: get your track list
Steps I have taken...

I first ran this command: 
```
python3 filter_tracklist.py track_list_2MflatPt_regression.txt 2.0 200000 track_list_200k_cotBeta_lt2.txt
```
After I ran the command I want to verify, so now I will do...
```
python3 check_tracklist_beta.py track_list_200k_cotBeta_lt2.txt

python plot_tracklist.py \
  track_list_200k_cotBeta_lt2.txt \
  --output track_list_200k_distributions.png
```
