# Installation

There are two ways to install *NiBabies*:
    - using container technologies; or
    - within a manually prepared environment, also known as *bare-metal*.

## Container Installation

Given its extensive dependencies, the easiest way to get up and running with *NiBabies* is by using a container service, such as [Docker](https://www.docker.com/get-started) or [Apptainer](https://apptainer.org/).

### Working with Docker

Images are hosted on our [Docker Hub](https://hub.docker.com/r/nipreps/nibabies).
To pull an image, the specific version tag must be specified in order to pull the images.
For example, to pull the first release in the 24.0.0 series, you can do:

```shell
docker pull nipreps/nibabies:24.0.0
```

There are also a few keyword tags, `latest` and `unstable`, that serve as special pointers.
`latest` points to the latest release (excluding any betas or release candidates).
`unstable` points to the most recent developmental change, and should only be used to test new features or fixes.

:::{tip}

Using versioned releases (such as `24.0.0`) are preferred to using the `latest` tag, as they are more explicit.
:::

### Working with Apptainer (formerly Singularity)

The easiest way to create an Apptainer image is to build from the [Docker](#working-with-docker) images hosted online.

```bash
apptainer build nibabies-version.sif docker://nipreps/nibabies:version
```

## Installing the nibabies-wrapper

The `nibabies-wrapper` is a lightweight Python tool to facilitate running `nibabies` within a container service.
To install or upgrade to the current release:
```
$ pip install --update nibabies-wrapper
```

For further details, see [](usage.md#using-the-nibabies-wrapper).

## Bare-metal Installation

If you would prefer to install this tool natively, you can refer the [Dockerfile](https://github.com/nipreps/nibabies/blob/master/Dockerfile) as a guide for all the dependencies.
