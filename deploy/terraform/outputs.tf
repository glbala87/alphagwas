# AlphaGWAS Terraform Outputs

output "vpc_id" {
  description = "VPC ID"
  value       = module.vpc.vpc_id
}

output "private_subnets" {
  description = "Private subnet IDs"
  value       = module.vpc.private_subnets
}

output "public_subnets" {
  description = "Public subnet IDs"
  value       = module.vpc.public_subnets
}

output "ecr_repository_url" {
  description = "ECR repository URL"
  value       = aws_ecr_repository.alphagwas.repository_url
}

output "ecs_cluster_name" {
  description = "ECS cluster name"
  value       = aws_ecs_cluster.main.name
}

output "ecs_cluster_arn" {
  description = "ECS cluster ARN"
  value       = aws_ecs_cluster.main.arn
}

output "alb_dns_name" {
  description = "ALB DNS name"
  value       = aws_lb.main.dns_name
}

output "alb_zone_id" {
  description = "ALB Route53 zone ID"
  value       = aws_lb.main.zone_id
}

output "api_url" {
  description = "API URL"
  value       = "http://${aws_lb.main.dns_name}"
}

output "s3_bucket_name" {
  description = "S3 data bucket name"
  value       = aws_s3_bucket.data.id
}

output "s3_bucket_arn" {
  description = "S3 data bucket ARN"
  value       = aws_s3_bucket.data.arn
}

output "cloudwatch_log_group" {
  description = "CloudWatch log group name"
  value       = aws_cloudwatch_log_group.alphagwas.name
}

output "rds_endpoint" {
  description = "RDS endpoint"
  value       = var.enable_rds ? aws_db_instance.main[0].endpoint : "N/A (RDS disabled)"
}

output "rds_database_name" {
  description = "RDS database name"
  value       = var.enable_rds ? aws_db_instance.main[0].db_name : "N/A (RDS disabled)"
}

# Deployment instructions
output "deployment_instructions" {
  description = "Instructions for deploying the application"
  value       = <<-EOT

    AlphaGWAS has been deployed to AWS!

    Next steps:

    1. Build and push Docker image:
       aws ecr get-login-password --region ${var.aws_region} | docker login --username AWS --password-stdin ${aws_ecr_repository.alphagwas.repository_url}
       docker build -t ${aws_ecr_repository.alphagwas.repository_url}:latest .
       docker push ${aws_ecr_repository.alphagwas.repository_url}:latest

    2. Update ECS service to deploy:
       aws ecs update-service --cluster ${aws_ecs_cluster.main.name} --service ${var.project_name}-api --force-new-deployment

    3. Access the API:
       curl http://${aws_lb.main.dns_name}/health

    4. View logs:
       aws logs tail /ecs/${var.project_name} --follow

  EOT
}
