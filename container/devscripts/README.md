# devscripts
development scripts for image maintenance

# notes

### building docker images
- user must be in `docker` group, `pwd` must contain `docker-compose.yml`
```
docker-compose build ci
docker-compose build dev
```

### pushing to dockerhub
example for CI image: give it the dockerhub repo name, with the tag `ci`
```
docker login
docker images
docker tag largexeic-ci:latest cjdilks/largex-eic:ci
docker push !$
docker logout
```
