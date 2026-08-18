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
## Step 2: PixelAV dataset generation
```
gcc -O2 ppixelav2_list_trkpy_n_2f_10ps_48x192.c -msse -lm -o ppixelav2_list_trkpy_n_2f_10ps_48x192
```
then after...

```
./ppixelav2_list_trkpy_n_2f_10ps_48x192 <frun> <runsize>
# frun    = run index (1–TEMPMAX); selects which block of track_list.txt to process
# runsize = tracks per run (we use 5000; max = NMUON = 50000)
# Output file: pixel_clusters_d<filebase+frun>.out  (filebase read from ppixel2.init)
```
