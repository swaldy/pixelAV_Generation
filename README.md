# pixelAV_Generation
```
2M original tracks
       │
       │ filter_tracklist.py
       │ |cotBeta| < 2
       │ select 200,000
       ▼
track_list_200k_cotBeta_lt2.txt
       │
       │
       ▼
PixelAV 48×192, 200 ps
       │
       ├── block 1:  tracks 1–5,000
       ├── block 2:  tracks 5,001–10,000
       ├── ...
       └── block 40: tracks 195,001–200,000
```       
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
