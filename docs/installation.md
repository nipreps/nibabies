# Installation

## Container Installation

Given its extensive dependencies, the easiest way to get up and running with *NiBabies* is by using a container service, such as [Docker](https://www.docker.com/get-started) or [Singularity](https://sylabs.io/singularity/).

### Working with Docker

Images are hosted on our [Docker Hub](https://hub.docker.com/r/nipreps/nibabies).
To pull an image, the specific version tag must be specified in order to pull the images.
For example, if you want to pull version `21.0.0`, you would use the following command.
```
$ docker pull nipreps/nibabies:21.0.0
```

There are also a few keyword tags, `latest` and `unstable`, that serve as special pointers.
`latest` points to the latest release (excluding any betas or release candidates).
`unstable` points to the most recent developmental change, and should only be used to test new features or fixes.

### Working with Singularity

The easiest way to create a Singularity image is to build from the [Docker](#Working-with-Docker) images hosted online.
For example, if you want to build version `21.0.0`, you would use the following command.
```
$ singularity build nibabies-21.0.0.sif docker://nipreps/nibabies:21.0.0
```

### Installing the nibabies-wrapper

The `nibabies-wrapper` is a lightweight Python tool to facilitate running `nibabies` within a container service.
To install the current release:
```
$ pip install nibabies-wrapper
```

You can find the [usage instructions here](./usage.md#Using-the-nibabies-wrapper)

## Bare-metal Installation

If you would prefer to install this tool natively, you can refer the [Dockerfile](./Dockerfile) as a guide for all the dependencies.