@echo off
echo removing old files..
rm -r docs/*

echo building new...
mkdocs build --config-file mkdocs/mkdocs.yml

echo copying files..
cp -r mkdocs/site/** docs

echo opening site...
"docs/index.html"
echo done!