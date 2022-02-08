This is the code for computing psuedo Cl from the masked/weighted shear in flat approximation.
The code includes the following files

 * README
 * Makefile
 * iparam.sh : a shell script to set parameters (grid number, box size...)
 * iparam.f90 : set basic parameters
 * default.inp : default parameters for iparam.f90
 * test.sh: a shell script to test the shear power reconstruction
 * sheardist.f90 : make a shear catalog from Gaussian shear fields
 * starmask.f90 : simulate mask of weak lensing data 
 * maskinfo.f90 : set parameters for simulated masks
 * modemat_square.f90 : compute mode coupling matrix of a square patch of sky
 * modemat_inside.f90 : compute mode coupling matrix of the simulated mask
 * shearpow.f90 : compute the shear power spectrum from shear catalogs
 * deconv.f90 : (de)convolve the shear power spectrum
 * collect.f90 : collect files of output power spectrum into one file
 * sub.f90 : subroutines

To run the code, 

 * "make setparam" to set parameters (e.g., grid number, box size, binnings)
 * "make all" to compile codes 
 * "make test" to test reconstruction method using Gaussian shear fields 
   it may take about 5mins
 
You can see the results of binned power spectrum in the file "sumpow.dat"
 * 1st line: l
 * 2nd line: input power (convergence)
 * 3rd line: deconvolved power (E-mode shear)
 * 4th line: convolved power (E-mode shear)
 * 5th line: masked power (E-mode shear)
 * 6-9th line: same as 2-5th line but for B-mode shear power

# 内容
本コードは2次元平面として近似できる天域領域において、弱い重力レンズ効果のゆがみ(シアー)場のパワースペクトルを測定するコードです。

弱い重力レンズ効果は宇宙大規模構造による重力の効果によって遠方銀河の形がゆがむ現象で、宇宙のダークマター分布を調べることができるユニークな観測量です。
銀河はもともと真円ではなく、弱い重力レンズ効果は楕円率をほんの数％程度変える程度の小さい効果であるため、
広い天域にわたって遠方にある多くの銀河の形を精密に観測し、統計的に重力レンズゆがみの大きさを引き出す処理が必要になります。
しかし、天域の方角によって精度良く観測できる銀河の数が少なかったり、明るい星が手前にあるとその背後にある銀河の形を測れないことがあり、
重力レンズゆがみの場は一様ではありません。

擬似スペクトル法は、場の非一様性の影響を計算し、本来の重力レンズゆがみの角度パワースペクトルを測定する方法です。
以下は、レイトレーシングシミュレーションから実際の観測を想定した非一様な重力レンズゆがみ場を作成し、
シミュレーションに入力した重力レンズパワースペクトル(黒の実線)が再現できることを確かめた図です。
非一様性の影響を全く考慮しないと、上パネルの青点のようにパワースペクトルの振幅が大きくずれます。
擬似スペクトルの手法を用いることで(上パネルの赤点)、正しいスペクトル(黒の実線)を再現することができます。
実際の銀河のように形がもともとゆがんでいる(平均的な楕円率の大きさが22％程度)場合でも、正しいスペクトルを再現することができます(右図)

また重力レンズによるゆがみのパターンはEモードとよばれる放射状の成分のみで、渦状のBモード成分は出てきませんが、
観測データの非一様性を無視して計算してしまうとBモードが出てきてしまいます。擬似スペクトル法を使うことで観測的なBモードの影響を
小さく抑えることができます(図の下パネル)

<img width="623" alt="スクリーンショット 2021-10-14 14 02 26" src="https://user-images.githubusercontent.com/86592645/137254998-644ab80a-8409-45a9-b4b3-9b6556aa7e7a.png">

本研究成果はすばるハイパーシュプリームカムの観測データに応用し、宇宙論情報を引き出すのにも使われました。

# References
- Shear Power Spectrum Reconstruction using Pseudo-Spectrum Method  
Chiaki Hikage, Masahiro Takada, Takashi Hamana, David Spergel  
Mon. Not. Roy. Astron. Soc., Vol.412 Issue 1 (2011), 65-74

- Cosmology from cosmic shear power spectra with Subaru Hyper Suprime-Cam first-year data  
Chiaki Hikage, Masamune Oguri, Takashi Hamana, et al.  
Publ. Astron. Soc. Japan, Vol.71, Issue 2, 43 (2019)  
