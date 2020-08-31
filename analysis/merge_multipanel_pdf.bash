cp pdf_density_weighted_isocool_tctf_0.3_cr_1.00_transport_pdf_compare.png img1.png
cp pdf_temperature_weighted_isocool_tctf_0.3_cr_1.00_transport_pdf_compare.png img2.png
cp pdf_creta_weighted_isocool_tctf_0.3_cr_1.00_transport_pdf_compare.png img3.png

convert -resize 98% img2.png img2.png
mogrify -background white -gravity North -extent 2337x3354 img2.png
convert -crop 96.5x0%x+0+0 -gravity West img2.png img2.png

convert -crop 96.5x0%x+0+0 -gravity East img3.png img3.png
convert -resize 99.4% img3.png img3.png
mogrify -background white -gravity West -extent 2467x3016 img3.png

convert +append -gravity South img1.png  img2.png img3.png merged.png

