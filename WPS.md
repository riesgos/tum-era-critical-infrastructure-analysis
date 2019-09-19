This is (the start of) a documentation about the WPS integration
================================================================

The system reliability service is included at the moment on the GFZ
RIESGOS WPS Repository, you can find here:

http://rz-vm140.gfz-potsdam.de/wps/WebProcessingService

The current techology used there is
1) the WPS Server by 52°North (https://github.com/52North/WPS)
2) the GFZ-command-line-repository (
https://github.com/riesgos/gfz-command-line-tool-repository)


Where to find all necessary elements?
-------------------------------------

To set it up yourself you can use the following docker image:
https://hub.docker.com/r/gfzriesgos/riesgos-wps

You also need docker image for the system reliability service:

https://hub.docker.com/r/gfzriesgos/system_reliability

The last thing you need is the json configuration file for the
service, that is provided here:

https://github.com/riesgos/gfz-command-line-tool-repository/blob/development-sprint/src/main/resources/org/n52/gfz/riesgos/configuration/system-reliability.json

(Please note that this is a branch for the development sprint we did
in Oberpfaffenhofen in August 2019; this branch may be merged to the
master later).


How to set up it yourself
-------------------------

- pull the docker image for the system reliability service
- pull the docker image for the gfz-command-line-repository
- start the gfz-command-line-repository container
- change the configuration (servername, passwords, ...)
- copy the json configuration file into the docker-container for the
  running gfz-command-line-repository in the folder
  /usr/share/riesgos/json-configurations/

Any problem on setting up the service?
--------------------------------------

You can find further information on how it works in the documentation
 https://github.com/riesgos/gfz-command-line-tool-repository
 
Maybe it is also helpful to read the commands in the dockerfile for
the gfz-command-line-repository, for example this one:
https://hub.docker.com/r/gfzriesgos/riesgos-wps/dockerfile

because this clears that there is the need to pass the socket for the
docker process to the server, so that the WPS server can run docker
commands.
 
It there is still something not clear, please contact
matthias.ruester@gfz-potsdam.de or nils.brinckmann@gfz-potsdam.de,
so that we can improve the documentation.

How does it works?
------------------

The server and the gfz-command-line-repository are running inside of
docker, containing a tomcat, a geoserver and the wps server.
On top of this is the jar with the command line repository which
provides a skeleton to run command line processes as wps services.

The services themselfes run also inside of docker to make them
independent from each other in terms of seperated runs (no conflicts
with intermediate temporary files) and dependencies.

Any kind of command line program can be executed, no matter if it is
a python script, a bash script, a c/c++ program or even fortran code.
Any combinations of those are possible too.

The processes themselves are configured via json files, that explains
how to call the process, how to provide input (writing in stdin,
giving command line arguments or writing files before the program
starts).
Same is true for the output, so the data can be read from stderr,
stdout or output files.

The command-line-repository can be configured where to search for
those json files and integrates them on runtime.

By default this folder is in the docker container under:
/usr/share/riesgos/json-configurations/

What does it mean for the system reliability service?
-----------------------------------------------------

We wrote the dockerfile with a base image and the dependencies for the
script. Those include python3, scipy, numpy and networkx.

The json configuration says that we use the
```
System_Reliability_Analysis_local.py
```
file to execute the overall program.

The script takes three input parameter:
- a file for the intensity values (in form of shakemaps)
- a name for the country to use (chile / ecuador)
- a name for the hazard to use for the fragility functions (earthquake
  / lahar)
  
(At the time of writing only the earthquake configuration works; there
is currently no implementation for reading the intensity / several
intensities from lahars).

We also set a default output name as the script uses different ones
depending on the chosen country.
