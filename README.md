# HLA-DPB1 Expression/TCE Suite

# Table of Contents
- [Background](#background-)
- [Suite](#suite-)
    - [DPdb](#DPdb-)
    - [DP microservice](#dp-microservice-)
    - [DP tool](#dp-tool-)
- [DP Tool Rest End Points Docker Container](#dp-tool-rest-end-points-docker-container-)
    - [Prerequisite](#prerequisite-)
    - [Docker image build](#docker-image-build-)
    - [Docker container launch](#docker-container-launch-)
    - [Stopping the running container](#stopping-the-running-container-)
- [DP Tool Web Application Docker Container](#dp-tool-web-application-docker-container-)
    - [Webapp Docker Prerequisite](#webapp-docker-prerequisite-)
    - [Webapp Docker image build](#webapp-docker-image-build-)
    - [Webapp Docker container launch](#webapp-docker-container-launch-)
    - [Stopping the running webapp container](#stopping-the-running-webapp-container-)
- [DP Tool Deployment With Docker Compose](#dp-tool-deployment-with-docker-compose-)
    - [Deployment Prerequisite](#deployment-prerequisite-)
    - [App Deployment](#app-deployment-)
    - [App Undeploy](#app-undeploy-)
- [DP Tool Production Deployment](#dp-tool-production-deployment-)
    - [Unified Container Deployment](#unified-container-deployment-)
    - [Decoupled Containers Deployment](#decoupled-containers-deployment-)


# Background [⤴](#table-of-contents)
This repository supports the annotation of *HLA-DPB1* from reference data on the [IPD-IMGT/HLA database](https://www.ebi.ac.uk/ipd/imgt/hla/) via two models as detailed by [Sajulga et al. 2023](https://pubmed.ncbi.nlm.nih.gov/37126658/):

| Model | Background | Reference data | Papers |
| -     | -          | -              | -      |
| [T-Cell Epitope (TCE)](https://www.ebi.ac.uk/ipd/imgt/hla/dpb.html) | Amino acid motifs on the antigen-recognition domain (encoded by exon 2) that categorize the 1,800+ alleles into different TCE groups (1 - most immunogenic > 2 > 3 - least immunogenic). | [IPD-IMGT/HLA](https://github.com/ANHIG/IMGTHLA/blob/Latest/tce/dpb_tce.csv) | [Crivello et al. 2015](https://www.astctjournal.org/article/S1083-8791(14)00656-9/fulltext)
| Expression level | Based off of the 3' UTR expression marker in high linkage disequilibrium with exon 3. Patients that have mismatched high expression alleles have higher risks for acute graft-versus-host disease. | [IPD-IMGT/HLA](https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/alignments/DPB1_nuc.txt) | [Petersdorf et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4560117/), [Schöne et al. 2018](https://dx.doi.org/10.1016%2Fj.humimm.2017.11.001), [Petersdorf et al. 2020](https://ascopubs.org/doi/full/10.1200/JCO.20.00265)


There are several methods of accessing this information and annotating *HLA-DPB1* alleles, genotypes, and recip-donor pairs.
| Name            | Type              | Application     |
| -               | -                 | -               |
| [DPdb](#DPdb-)            | Python package    | Incorporation into other algorithms or scripts |
| [DP microservice](#dp-microservice-) | REST microservice | Same as DPdb, but easier for integration with other programming languages or platforms. |
| [DP tool](#dp-tool-)         | Angular web interface | No programming needed. For single-use, easy annotation via browser |

# Suite [⤴](#table-of-contents)
## DPdb [⤴](#table-of-contents)

### Setup

To begin, ensure that you have Python3 installed. To check, issue this command to verify your python version:
```
python --version
```

If Python3 is not installed, please download it from [here](https://www.python.org/downloads/).

#### Virtual environment

If Python3 is readily available, set up your virtual environment by running these commands:
```
python3 -m venv venv
source venv/bin/activate
```

#### Python package installation
```
pip install hlann
```

Pip is the package installer for Python. It comes pre-packaged with Python. 


#### Python package loading
Open a python instance.
```
python
```

Load the python package for annotation.
```
from hlann.hlann import HLAnn
dp_db = HLAnn(db_version='3540', verbose=True)
```

### Usage
#### Allele
Command
```
dp_db.annotate_allele('DPB1*01:AETTG')
```
<details>
<summary>Result</summary>
<pre>
{'name': 'DPB1*01:AETTG',
 'tce': '3',
 'resolution': 'intermediate',
 'matched': False,
 'expr_level': '~high',
 'expr_annot_type': '~experimental',
 'alleles': [{'allele_expr': 'DPB1*01:01',
   'CIWD_TOTAL': 'C',
   'exon3_motif_ref': 'ACCACTC',
   'utr3_motif_ref': 'G',
   'tce': '3',
   'resolution': 'high',
   'expr_level': 'high',
   'experimental': True},
  {'allele_expr': 'DPB1*162:01',
   'CIWD_TOTAL': 'WD',
   'exon3_motif_ref': 'GTTGTCT',
   'utr3_motif_ref': 'A',
   'tce': '3',
   'resolution': 'high',
   'expr_level': 'low',
   'experimental': False},
  {'allele_expr': 'DPB1*417:01',
   'CIWD_TOTAL': 'WD',
   'exon3_motif_ref': 'ACCACTC',
   'utr3_motif_ref': 'unknown',
   'tce': '3',
   'resolution': 'high',
   'expr_level': 'high',
   'experimental': False}]}
</pre>
</details>
<br>

#### Genotype
Command
```
dp_db.annotate_genotype('DPB1*01:AETTA+DPB1*04:AETTB')
```
<details>
<summary>Result</summary>
<pre>
{
  "allele_one": {
    "alleles": [
      {
        "CIWD": "C",
        "allele_expr": "DPB1*01:01",
        "exon3_motif_ref": "ACCACTC",
        "experimental": true,
        "expr_level": "high",
        "resolution": "high",
        "tce": "3",
        "utr3_motif_ref": "G"
      },
      {
        "CIWD": "WD",
        "allele_expr": "DPB1*417:01",
        "exon3_motif_ref": "ACCACTC",
        "experimental": false,
        "expr_level": "high",
        "resolution": "high",
        "tce": "3",
        "utr3_motif_ref": "unknown"
      }
    ],
    "expr_annot_type": "experimental",
    "expr_level": "high",
    "matched": false,
    "name": "DPB1*01:AETTA",
    "resolution": "intermediate",
    "tce": "3"
  },
  "allele_two": {
    "alleles": [
      {
        "CIWD": "C",
        "allele_expr": "DPB1*04:01",
        "exon3_motif_ref": "GTTGTCT/GTTGCCT/GTTATCT",
        "experimental": true,
        "expr_level": "low",
        "resolution": "high",
        "tce": "3/0",
        "utr3_motif_ref": "A"
      },
      {
        "CIWD": "C",
        "allele_expr": "DPB1*126:01",
        "exon3_motif_ref": "GTTGTCT",
        "experimental": false,
        "expr_level": "low",
        "resolution": "high",
        "tce": "3",
        "utr3_motif_ref": "A"
      },
      {
        "CIWD": "I",
        "allele_expr": "DPB1*350:01",
        "exon3_motif_ref": "ACCACTC",
        "experimental": false,
        "expr_level": "high",
        "resolution": "allelic",
        "tce": "3",
        "utr3_motif_ref": "unknown"
      },
      {
        "CIWD": "WD",
        "allele_expr": "DPB1*415:01",
        "exon3_motif_ref": "GTTGTCT",
        "experimental": false,
        "expr_level": "low",
        "resolution": "allelic",
        "tce": "3",
        "utr3_motif_ref": "A"
      }
    ],
    "expr_annot_type": "~experimental",
    "expr_level": "~low",
    "matched": false,
    "name": "DPB1*04:AETTB",
    "resolution": "intermediate",
    "tce": "3"
  },
  "genotype": "DPB1*01:AETTA+DPB1*04:AETTB"
}
</pre>
</details>
<br>


#### Matches
Command
```
dp_db.annotate_match('DPB1*04:01+DPB1*40:01', 'DPB1*40:01+DPB1*40:01')
```
<details>
<summary>Result</summary>
<pre>
{
  "directionality": "GvH",
  "genotype_donor": {
    "allele_one": {
      "alleles": null,
      "expr_annot_type": "experimental",
      "expr_level": "low",
      "matched": true,
      "name": "DPB1*40:01",
      "resolution": "high",
      "tce": "3"
    },
    "allele_two": {
      "alleles": null,
      "expr_annot_type": "experimental",
      "expr_level": "low",
      "matched": true,
      "name": "DPB1*40:01",
      "resolution": "high",
      "tce": "3"
    },
    "genotype": "DPB1*40:01+DPB1*40:01"
  },
  "genotype_recipient": {
    "allele_one": {
      "alleles": null,
      "expr_annot_type": "experimental",
      "expr_level": "low",
      "matched": false,
      "name": "DPB1*04:01",
      "resolution": "high",
      "tce": "3"
    },
    "allele_two": {
      "alleles": null,
      "expr_annot_type": "experimental",
      "expr_level": "low",
      "matched": true,
      "name": "DPB1*40:01",
      "resolution": "high",
      "tce": "3"
    },
    "genotype": "DPB1*04:01+DPB1*40:01"
  },
  "grade": "MA",
  "matched_alleles_don": [
    "DPB1*40:01",
    "DPB1*40:01"
  ],
  "matched_alleles_pat": [
    "DPB1*40:01"
  ],
  "mismatched_allele_pat_expr_level": "low",
  "mismatched_alleles_don": [],
  "mismatched_alleles_pat": [
    "DPB1*04:01"
  ],
  "tce_match": "Permissive"
}
</pre>
</details>
<br>



## DP microservice [⤴](#table-of-contents)
### Bootstrapping

To begin, ensure that you have Python3 installed. To check, issue this command to verify your python version:
```
python --version
```

If Python3 is not installed, please download it from [here](https://www.python.org/downloads/).

If Python3 is readily available, set up your virtual environment by running these commands:
```
python3 -m venv venv
source venv/bin/activate
```

Pip is the package installer for Python. It comes pre-packaged with Python. This will be used to install our requirements as such:
```
pip install --upgrade pip
pip install -r requirements.txt
```

Once installed, *behave* will be available for testing and *Flask* will be available to set up the web service.

### Initialization

Initialize the web service via this command:
```
python server.py
```

Once initialized, you may use the REST API endpoints.

> *REST → **Re**presentational **S**tate **T**ransfer*
 *API → **A**pplication **P**rogramming **I**nterface*

For example, you can retrieve data from an endpoint by using cURL as shown below.


### Testing

Running all the tests in this repository is as simple as running this command:
```
behave
```

## DP tool [⤴](#table-of-contents)

This front-end graphical user interface (GUI) was created using Angular 8.

### Bootstrapping

We will need to go into the web app project's root folder
```
cd webapp
```

Since our web application uses JavaScript (Angular 8), install Node.js (≥10.9) and npm (node package manager) [here](https://nodejs.org/en/download/) if ```npm``` is not a recognized command in your terminal.

Through npm, we can install our dependencies by running:
```
npm install
```

You can also install the Angular CLI (to run `ng` `<command>` as seen in the next section) via:
```
npm install -g @angular/cli
```

### Initialization

Once finished, ensure that the back-end REST server has been initialized on http://0.0.0.0:5010/ as detailed [here](#dp-microservice-).

And then run a local development server:
```
ng serve
```

The web application will now be available on https://0.0.0.0:4200/.

## DP Tool Rest End Points Docker Container [⤴](#table-of-contents)

### Prerequisite [⤴](#table-of-contents)
The containerization is facilitated by [Docker Container](https://www.docker.com/resources/what-container). 

To be able to run Docker container, a docker set up and configuration is necessary. The installation details can be found in the [official docker documentation](https://docs.docker.com/get-docker/).

### Docker image build [⤴](#table-of-contents)
To build the image, navigate to the directory where `Dockerfile-flask` is located
Execute the comand (keep an eye on required "." at the end of the command)
```
docker build -t nmdpbioinformatics/dp-tool-backend:latest -f Dockerfile-flask .
```
Now the image should be built and available in the local docker registry

### Docker container launch [⤴](#table-of-contents)
To start a container form docker image (built in the last step) we need to execute the following command

```
docker run -d -p 5010:5010 nmdpbioinformatics/dp-tool-backend:latest
```
Upon successful execution a container id should comeout. We can see the container if it is up by executing
```
docker ps -a
```
That should show us if the container is up and running. If it is up then we should be able to navigate to http://localhost:5010/ to see the API landing page.

### Stopping the running container [⤴](#table-of-contents)
We have to obtain the container id by executing 
```
docker ps -a
```
Then we have to execute 
```
docker stop $CONTAINER_ID
```

## DP Tool Web Application Docker Container [⤴](#table-of-contents)

### Webapp Docker Prerequisite [⤴](#table-of-contents)

[DP Tool Rest End Points Docker Container](#dp-tool-rest-end-points-docker-container-) is up and runninng for the webapp container to be working properly (for build that is not necessary).

The containerization is facilitated by [Docker Container](https://www.docker.com/resources/what-container). 

To be able to run Docker container, a docker set up and configuration is necessary. The installation details can be found in the [official docker documentation](https://docs.docker.com/get-docker/).

### Webapp Docker image build [⤴](#table-of-contents)
To build the image, navigate to the `webapp` directory where `Dockerfile` is located.
We will need to go into the web app project's root folder using the following command
```
cd webapp
```
Execute the command below(keep an eye on required "." at the end of the command)
```
docker build --build-arg CONFIGURATION="" -t nmdpbioinformatics/dp-tool-ui-app .
```
Now the image should be built and available in the local docker registry.

### Webapp Docker container launch [⤴](#table-of-contents)
To start the webapp container form docker image (built in the last step) we need to execute the following command

```
docker run -d -p 80:80 -t nmdpbioinformatics/dp-tool-ui-app:latest
```
Upon successful execution a container id should comeout. We can see the container if it is up by executing by executing
```
docker ps -a
```
That should show us if the container is up and running. If it is up then we should be able to navigate to http://localhost:80/ to see the Web App landing page.

### Stopping the running webapp container [⤴](#table-of-contents)
We have to obtain the container id by executing 
```
docker ps -a
```
Then we have to execute 
```
docker stop $CONTAINER_ID
```

## DP Tool Deployment With Docker Compose [⤴](#table-of-contents)
**Warning: This segment is designed for local development purpose, seperate `docker-compose` file configuration is needed be developed for production deployment with docker-compose.**
The app front end and back end both can be deployed using [Docker Compose](https://docs.docker.com/compose/).

### Deployment Prerequisite [⤴](#table-of-contents)
[Docker Compose](https://docs.docker.com/compose/) have to be installed and docker registry have to contain the images a) `be-the-match/dp-tool-ui-app` and b) `be-the-match/dp-tool-backend`.

### App Deployment [⤴](#table-of-contents)
Simply execute the following command where `docker-compose.yml` file is located to deploy the application frontend and backend
```
docker-compose up -d
```

### App Undeploy [⤴](#table-of-contents)
To undeploy the app simply execute the following command where `docker-compose.yml` file is located
```
docker-compose down
```

## DP Tool Production Deployment [⤴](#table-of-contents)
The production deployment has two model, a) unified container model and b) segregated contaienrs model. 

a) The unified container packes both the backend (python-flask-gunicorn) and front end (angular and nginx) into one docker docker images while during the runtime Nginx acts as an webserver for front end and reverse proxy for backend. 

b) The segregated container deployment would provide independent scalling of backend and front end cluster should there be any need for it. Although the decoupling might be desireable under certain circumstances but this feature would require setting up an approprite network using docker compose or Kubernetes and are currently not available.

### Unified Container Deployment [⤴](#table-of-contents)
The production `apiUrl` should be adjusted with correct server in the file `webapp/src/environments/environment.prod.ts`

To build the docker image the following command may be executed in the project root directory:
```
docker build --build-arg CONFIGURATION="production" -t nmdpbioinformatics/dp-tool .
```
After successful build we should have the docker image available in our local docker registry.

To deploy the app now we can use the following command. The application should be available in your domain, i.e. `http://host:80/`
```
docker run -d -p 80:80 -t nmdpbioinformatics/dp-tool:latest
```
To stop the app container, We have to obtain the container id by executing 
```
docker ps -a
```
Then we have to execute 
```
docker stop $CONTAINER_ID
```

### Decoupled Containers Deployment [⤴](#table-of-contents)
This feature is under development now.