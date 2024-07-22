cp ~/projects/.fossils/pipeline.fossil .
podman build . --tag buildpipeline
echo podman run --name pipeline_builder --replace -v ./dist:/out localhost/buildpipeline:latest bash scripts/build_python3.9.sh
podman run --name pipeline_builder --replace -v ./dist:/out localhost/buildpipeline:latest bash scripts/build_python3.10.sh
echo podman run --name pipeline_builder --replace -v ./dist:/out localhost/buildpipeline:latest bash scripts/build_python3.11.sh
echo podman run --name pipeline_builder --replace -v ./dist:/out localhost/buildpipeline:latest bash scripts/build_python3.12.sh
rsync dist/* fermat:pkgs/
