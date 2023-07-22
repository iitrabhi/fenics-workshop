## Prerequisites

To follow along with the examples you need to install docker on your system. You need Windows 10/11 Education or Professional for this to work. This does not work on Windows 10/11 Home.

* [Docker](https://www.docker.com/products/docker-desktop)
* [Paraview](https://www.paraview.org/download/)
* [CMDER](https://cmder.net/) (Only for Windows)

After installation open `cmder` and then go-to Settings(Win+Alt+P)➡import and choose the `cmlab.xml` provided in the repository.

## Install FEniCS

Once the docker system in installed and running open CMDER/terminal and run:

```
docker pull iitrabhi/fenics_notebook
```

## Running

To start the notebook server use:

```
docker run -p 8888:8888 -v host_system_path:/root/ -w /root/ iitrabhi/fenics_notebook
```

Note: you should replace the variable `host_system_path` with the path of the folder that contains your code. e.g. If  `D:\Codes` contains your code then to start the command line interface you have to run:

```
docker run -p 8888:8888 -v D:\Codes:/root/ -w /root/ iitrabhi/fenics_notebook
```

Once you run the above command in `cmder` you will get a URL starting with `http://127.0.0.1:8888/lab`. Press Control and click on the URL to open a new Jupyter notebook.

## If you have Windows home

https://fem-on-colab.github.io/packages.html