# devscripts
development scripts for image maintenance

# notes

### building docker images
- user must be in `docker` group, `pwd` must contain `docker-compose.yml`
```
docker-compose build dev
container/devscripts/dockerShell.sh   # for testing
```

### pushing to dockerhub
example for CI image: give it the dockerhub repo name, with the tag `dev`
```
docker login
docker images
docker tag sidiseic_dev:latest cjdilks/sidis-eic:latest
docker scan cjdilks/sidis-eic:latest  2>&1 | tee vulnerabilities.txt
docker push cjdilks/sidis-eic:latest
docker logout
```
