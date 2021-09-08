set -e

echo "Building mtRNA finder..."

cp -r $SRC_DIR/* $PREFIX/

mkdir $PREFIX/bin
cd $PREFIX/bin
ln -s  $PREFIX/mt_finder.py ./mt_finder
chmod +x ./mt_finder
