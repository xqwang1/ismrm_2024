# BART scripts to reconstruct wave-encoded multi-contrast / multi-echo data

# based on the BART tests function

# Authors: Xiaoqing Wang, Berkin Bilgic, Jose P. Marques, 2023-2024

# single-echo
bart wavepsf -x 640 -y 128 wave_psf
bart fft -iu 7 shepplogan_coil_ksp.ra img
bart resize -c 0 640 img wave_zpad
bart fft -u 1 wave_zpad wave_hyb
bart fmac wave_hyb wave_psf wave_acq
bart fft -u 6 wave_acq wave_ksp
bart wave coils.ra wave_psf wave_ksp reco




# multi-echo / contrast
bart wavepsf -x 640 -y 128 wave_psf
bart fft -iu 7 shepplogan_coil_ksp.ra img
bart repmat 5 6 img img_6echos 
bart resize -c 0 640 img_6echos wave_zpad
bart fft -u 1 wave_zpad wave_hyb_6echos
bart fmac wave_hyb_6echos wave_psf wave_acq
bart fft -u 6 wave_acq wave_ksp_6echos
bart wave coils.ra wave_psf wave_ksp_6echos reco_6echos



# LLR in the echo dimension 
bart wavepsf -x 640 -y 128 wave_psf
bart fft -iu 7 shepplogan_coil_ksp.ra img
bart repmat 5 6 img img_6echos 
bart resize -c 0 640 img_6echos wave_zpad
bart fft -u 1 wave_zpad wave_hyb_6echos
bart fmac wave_hyb_6echos wave_psf wave_acq
bart fft -u 6 wave_acq wave_ksp_6echos
bart wave -l -r1e-5 -b8 coils.ra wave_psf wave_ksp_6echos reco_6echos_llr
