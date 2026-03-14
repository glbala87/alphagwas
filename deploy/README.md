# AlphaGWAS Cloud Deployment

Deploy AlphaGWAS to AWS using Terraform.

## Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                         AWS Cloud                           │
│  ┌───────────────────────────────────────────────────────┐  │
│  │                       VPC                              │  │
│  │  ┌─────────────┐                                      │  │
│  │  │   Public    │  ┌─────────────┐                     │  │
│  │  │   Subnet    │  │     ALB     │◄── Internet         │  │
│  │  └─────────────┘  └──────┬──────┘                     │  │
│  │                          │                             │  │
│  │  ┌─────────────────────────────────────────────────┐  │  │
│  │  │               Private Subnet                     │  │  │
│  │  │  ┌──────────┐  ┌──────────┐  ┌──────────┐       │  │  │
│  │  │  │ ECS API  │  │ ECS API  │  │ ECS      │       │  │  │
│  │  │  │ Task 1   │  │ Task 2   │  │ Worker   │       │  │  │
│  │  │  └──────────┘  └──────────┘  └──────────┘       │  │  │
│  │  │                                                  │  │  │
│  │  │  ┌──────────────────┐  ┌───────────────────┐    │  │  │
│  │  │  │  RDS PostgreSQL  │  │   S3 Data Bucket  │    │  │  │
│  │  │  │   (optional)     │  │                   │    │  │  │
│  │  │  └──────────────────┘  └───────────────────┘    │  │  │
│  │  └─────────────────────────────────────────────────┘  │  │
│  └───────────────────────────────────────────────────────┘  │
└─────────────────────────────────────────────────────────────┘
```

## Prerequisites

- [Terraform](https://www.terraform.io/downloads.html) >= 1.0
- [AWS CLI](https://aws.amazon.com/cli/) configured
- [Docker](https://www.docker.com/) for building images

## Quick Start

### 1. Configure Variables

```bash
cd deploy/terraform
cp terraform.tfvars.example terraform.tfvars
# Edit terraform.tfvars with your settings
```

### 2. Initialize Terraform

```bash
terraform init
```

### 3. Plan Deployment

```bash
terraform plan
```

### 4. Deploy Infrastructure

```bash
terraform apply
```

### 5. Build and Push Docker Image

```bash
# Get ECR login
aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin <ECR_URL>

# Build and push
docker build -t <ECR_URL>:latest ..
docker push <ECR_URL>:latest
```

### 6. Deploy Application

```bash
# Force new deployment
aws ecs update-service --cluster alphagwas-cluster --service alphagwas-api --force-new-deployment
```

## Configuration Options

| Variable | Description | Default |
|----------|-------------|---------|
| `aws_region` | AWS region | `us-east-1` |
| `environment` | Environment (development/staging/production) | `development` |
| `api_cpu` | API CPU units | `512` |
| `api_memory` | API memory (MB) | `1024` |
| `api_desired_count` | Number of API containers | `2` |
| `enable_autoscaling` | Enable auto-scaling | `true` |
| `enable_rds` | Enable PostgreSQL database | `false` |

## Environments

### Development
- Single NAT gateway
- No deletion protection
- Minimal resources

### Production
- Multi-AZ NAT gateways
- Deletion protection enabled
- Auto-scaling enabled
- RDS with backups

```bash
# Production deployment
terraform apply -var="environment=production" -var="enable_rds=true"
```

## Monitoring

### CloudWatch Logs
```bash
aws logs tail /ecs/alphagwas --follow
```

### ECS Service Status
```bash
aws ecs describe-services --cluster alphagwas-cluster --services alphagwas-api
```

## Cleanup

```bash
terraform destroy
```

## Cost Estimates

| Environment | Monthly Cost (approx) |
|-------------|----------------------|
| Development | $50-100 |
| Production | $200-500 |

*Costs vary based on usage and data transfer.*

## Security

- All data encrypted at rest (S3, RDS)
- VPC with private subnets
- Security groups restrict access
- IAM roles with least privilege

## Troubleshooting

### Container not starting
```bash
aws logs tail /ecs/alphagwas --filter-pattern "ERROR"
```

### Health check failing
```bash
curl http://<ALB_DNS>/health
```

### Database connection issues
- Check security group rules
- Verify RDS endpoint in ECS task environment
