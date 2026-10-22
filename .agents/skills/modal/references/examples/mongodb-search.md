# Using MongoDB Atlas Vector and GeoJSON Search with Modal

This [example repo](https://github.com/modal-labs/search-california)
demonstrates how to use Modal and MongoDB together
to build a full-stack application.

The application is a hybrid search engine,
like the retrieval engines that power RAG chatbots,
but for satellite images of the state of California.
Images can be searched based on their
geospatial and temporal metadata or based on their semantic content
as captured by a pre-trained embedding model.

We use the [Clay foundation model](https://clay-foundation.github.io/model/index.html)
for embeddings and we source the images from the European Space Agency's
[Sentinel satellites](https://www.esa.int/Applications/Observing_the_Earth/Copernicus/The_Sentinel_missions).

You can take our deployment of the application for a spin
[here](https://modal-labs-examples--clay-hybrid-search.modal.run/).

## Overview

At the center of the application is a MongoDB Atlas instance
that stores metadata for a collection of satellite images.

Modal orchestrates the compute around that database:
retrieving data from elsewhere and storing it in the database,
computing vector embeddings for the data in the database,
and serving both a frontend and a client.

The dataflow looks something like this:

1. Every few days, the European Space Agency's
   [Sentinel Satellites](https://www.esa.int/Applications/Observing_the_Earth/Copernicus/The_Sentinel_missions)
   complete a full pass over the entire Earth, including California.
   The images are made available via a [public STAC API](https://element84.com/geospatial/introducing-earth-search-v1-new-datasets-now-available/).
2. Every day, we run a job on Modal that queries that STAC API
   for new images of California and store the metadata in a MongoDB Atlas
   database instance.
3. Asynchronously, we run a job on Modal to check which entries
   in the database don't have an associated embedding.
   These images are then sent to a serverless embedding service
   running on Modal. We send the resulting embeddings to the database.
4. We host a database client on Modal that allows the application's
   developers to manipulate the data. This client is also used by two
   Web Functions for vector and geospatial search queries powered by
   Atlas Search.
5. Finally, we run a simple static FastAPI server on Modal that serves
   an Alpine JS frontend for executing those queries and rendering their results.

This entire application —
from API queries and frontend UI to GPU inference and hybrid search —
is delivered using nothing but Modal and MongoDB Atlas.
Setting it up for yourself requires only credentials on these platforms
and a few commands, detailed below.

## Deploying the Backend

### Setup: Modal and MongoDB Atlas

You'll need a Python environment on your local machine.
Any recent version of Python should do.
Most of the dependencies will be installed in environments on Modal,
so you don't need to worry quite so much.

Follow the instructions [here](https://modal.com/docs/guide#getting-started)
to set up your Modal account.
The $30/month of compute included in Modal's free tier is
more than enough to deploy and host this example.

You'll also need an account on MongoDB Atlas.
You can find instructions [here](https://www.mongodb.com/docs/atlas/getting-started/).
We prefer the UI, rather than the CLI, for setup.
The free tier is more than sufficient to run this example.

You'll want to create a database called `modal-examples`.
Make sure it's accessible from [all IP addresses](https://stackoverflow.com/questions/66035947/allow-access-from-anywhere-mongodb-atlas).
In the process, you will create a database user with a password.
Navigate to the Modal Secrets dashboard [here](https://modal.com/secrets)
and add this information, as well as the connection string for your database,
to a Modal Secret based on the MongoDB template available in the dashboard.

### MongoDB Client (`database.py`)

If your Modal Secret and MongoDB Atlas instance are set up correctly,
you should be able to run the following command:

```bash
modal run -m backend.database::MongoClient.ping
```

Once that command is working, you can start manipulating the database
from Modal.

To start, you'll want to add an Area of Interest (AOI) to the database:

```bash
modal run -m backend.database --action add_aoi
```

By default, it's the state of California as defined by the GeoJSON
in this repository's `data` folder (originally retrieved from
[the `geojsonio` GitHub repository](https://github.com/ropensci/geojsonio/blob/7e4cc683ed3d6eec38a8cae5ce03fa6d82acafc7/inst/examples/california.geojson)).
You can pass a different GeoJSON file to the `add_aoi` action
with the `--target` flag.

The `modal run` command is used for one-off tasks.
To deploy the database client for use in other parts of the app
along with the webhooks that anyone can use to run search queries,
we use `modal deploy`:

```bash
modal deploy -m backend.database
```

Those webhooks come with interactive OpenAPI docs,
which you can access by navigating to the `/docs` route of the deployment's URL.
You should see that URL in the terminal output.
You can also find the URL in the app's [Modal dashboard](https://modal.com/apps).

For our deployment, the URL for the interactive docs for the geographic
search endpoint is
[`https://modal-labs-examples--clay-mongo-client-geo-search.modal.run/docs`](https://modal-labs-examples--clay-mongo-client-geo-search.modal.run/docs).

If you haven't yet run the backfill jobs for your database instance,
as described below, this search will not return any results,
but you can use it to check that the database client is deployed.

### Backfill and Updates (`extract.py`)

We add data to the database by querying the Sentinel STAC API for images.

Run the following command to search for images in the AOI
from the preceding week and add them to the database:

```bash
modal run -m backend.extract
```

You can either check the results via the Atlas UI
or by executing a search query in the database client's geo search webhook,
as described above.

To regularly update the database with new images,
we deploy the app defined in `extract.py`:

```bash
modal deploy -m backend.extract
```

This app also runs a regular job to add embeddings to the images
in the database.

But it doesn't compute the embeddings itself --
embeddings are provided by a separate service,
which is described next.

### Clay Embeddings Service (`embeddings.py`)

To build the environment for the embeddings service
and to test the embedding engine on some sample data,
execute the following command:

```bash
modal run -m backend.embeddings
```

To deploy this on Modal, we again use `modal deploy`:

```bash
modal deploy -m backend.embeddings
```

### Putting It All Together

Now that the embedding service is deployed,
we can add vectors by invoking the `enrich_vectors`
function in `extract` with `modal run`:

```bash
modal run -m backend.extract::enrich_vectors
```

This command will ensure all the images in the database have embeddings.

You should be able to observe them on records viewed via the Atlas UI
or by executing a search query via the database client's geo search webhook,
as described previously.

To use the embeddings for search, we recommend running the frontend UI,
which we walk through next.

## Deploying the Frontend

The frontend is much simpler than the backend.
It comprises a small Alpine JS app and a FastAPI Python server
to deliver it to client browsers.

You can play with our deployment of the frontend
[here](https://modal-labs-examples--clay-hybrid-search.modal.run/).

### Alpine App (`app.js`)

The Alpine app provides a basic interface for constructing geo search queries
by clicking on a map and viewing results.
Clicking on the returned images triggers a vector search for similar images.
Images can be furthermore filtered by date using the date pickers.

### FastAPI Server (`serve.py`)

This app is served to the client by a FastAPI server.

To deploy it, run the following command:

```bash
modal deploy -m frontend
```
