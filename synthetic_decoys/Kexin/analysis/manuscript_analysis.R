# manuscript analysis
# Neural net and ERT accuracy analysis
acc_NN = load_model_accuracy("Kexin/NN/nn_loo_acc_3.0.txt")
pdf(file = "Kexin_NN/nn_acc_3.0.pdf", width = 8.44, height = 6.96)
plot_color_bar(acc_NN, "acc")
dev.off()
pdf(file = "Kexin_NN/nn_recall1_3.0.pdf", width = 8.44, height = 6.96)
plot_color_bar(acc_NN, "recall1")
dev.off()

acc_ert = load_model_accuracy("Kexin/ERT/ert_loo_acc_3.0.txt")
pdf(file = "Kexin_RF/ert_acc_3.0.pdf", width = 8.44, height = 6.96)
plot_color_bar(acc_ert, "acc")
dev.off()
pdf(file = "Kexin_RF/ert_recall1_3.0.pdf", width = 8.44, height = 6.96)
plot_color_bar(acc_ert, "recall1")
dev.off()
