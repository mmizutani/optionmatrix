/* optionmatrix:                                                            */
/*                                                                          */
/* Options & Futures Matrix Modeler                                         */
/* View and Control Theoretical Option Chains                               */
/*                                                                          */
/* File: gtk_menu.cpp of optionmatrix                                       */
/*                                                                          */
/* Copyright (c) Anthony Bradford. 2012.                                    */
/* http://opensourcefinancialmodels.com                                     */
/* info@opensourcefinancialmodels.com                                       */
/*                                                                          */
/* optionmatrix may be freely redistributed.                                */
/* See file COPYING included with this distribution for license information */

/* 
   optionmatrix is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   optionmatrix program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <gtk/gtk.h>

#include "gtk_main.h"

#ifdef WIN32
#include <stdarg.h>
#include <windef.h>
#include <wingdi.h>
#include <winuser.h>
#include <Shellapi.h>
#endif

void on_menu_quit_activate(GtkWidget *widget, const struct _properties *properties)
{
  //g_print ("on_menu_quit_activate():\n");
  gtk_main_quit();
}

void on_menu_website_activate(GtkWidget *widget, const struct _properties *properties)
{

#ifndef WIN32
  GError *error = NULL;

  gtk_show_uri(NULL,PACKAGE_URL, GDK_CURRENT_TIME, &error);

  if (error != NULL)
  {
    if (error->domain == GDK_PIXBUF_ERROR)
      g_print ("GdkPixbufError: %s\n", error->message);
    else if (error->domain == G_FILE_ERROR)
      g_print ("GFileError: %s\n", error->message);
    else
      g_print ("An error \"%s\" in the unexpected domain:%d has occurred!\n", error->message, error->domain);

    g_error_free (error);
  }
#endif

#ifdef WIN32
  int ret;
  ret = (int) ShellExecute( NULL, "open", PACKAGE_URL, NULL, NULL, SW_SHOWNORMAL );
#endif

}

void on_menu_feedback_activate(GtkWidget *widget, const struct _properties *properties)
{

  g_print("on_menu_feedback_activate():\n");

#ifndef WIN32
  GError *error = NULL;

  gtk_show_uri(NULL,"mailto:info@opensourcefinancialmodels.com", GDK_CURRENT_TIME, &error);

  if (error != NULL)
  {
    if (error->domain == GDK_PIXBUF_ERROR)
      g_print ("GdkPixbufError: %s\n", error->message);
    else if (error->domain == G_FILE_ERROR)
      g_print ("GFileError: %s\n", error->message);
    else
      g_print ("An error \"%s\" in the unexpected domain:%d has occurred!\n", error->message, error->domain);

    g_error_free (error);
  }
#endif

#ifdef WIN32
  int ret;
  ret = (int) ShellExecute( NULL, "open", "mailto:info@opensourcefinancialmodels.com", NULL, NULL, SW_SHOWNORMAL );
#endif

}

static gboolean pngTimer(struct _properties *properties)
{
  g_print("pngTimer():\n");

  int index = 0;
  char imageName[1000] = { 0 };

  for( index = 0; index < 30; index++ )
  {
    sprintf(imageName, "images/%d.png", index + 1);
    g_print("Checking for name: %s\n", imageName);

    if( access( imageName, F_OK ) == -1 )
    {
      break;
    }
  }

  GdkPixbuf *logo;
  GError *error = NULL;

  if( index )
  {
    if( index == 1 )
    {
      logo = gdk_pixbuf_new_from_file("images/1.png", &error);
    } else
    {
      sprintf(imageName,"images/%d.png",1 + (int)(index * (rand() / (RAND_MAX + 1.0))));
      g_print("imageName = %s\n",imageName);
      logo = gdk_pixbuf_new_from_file (imageName, &error);
    }

    /* Set the application logo or handle the error. */
    if (error == NULL)
    {
      gtk_about_dialog_set_logo(GTK_ABOUT_DIALOG (properties->GtkInfo.dialogAbout), logo);

      g_object_unref(logo);

    } else
    {
        if (error->domain == GDK_PIXBUF_ERROR)
          g_print("GdkPixbufError: %s\n", error->message);
        else if (error->domain == G_FILE_ERROR)
          g_print("GFileError: %s\n", error->message);
        else
          g_print("An error in the unexpected domain:%d has occurred!\n", error->domain);

        g_error_free (error);
    }
  } else
  {
    g_print("No images available\n");
  }

  return GTK_IS_ABOUT_DIALOG(properties->GtkInfo.dialogAbout);
}

void on_menu_about_activate( GtkWidget *widget, struct _properties *properties )
{
  //g_print("on_menu_about_activate():\n");

  int index = 0;
  char imageName[1000] = { 0 };

  for( index = 0; index < 10; index++ )
  {
    sprintf( imageName, "images/%d.png", index + 1 );
    g_print("Checking for name: %s\n", imageName);

    if( access( imageName, F_OK ) == -1 )
    {
      break;
    }
  }

  const gchar *authors[] = {
    "Anthony Bradford",
    PACKAGE_BUGREPORT,
    PACKAGE_URL,
    "\nwith models from:",

#ifdef METAOPTIONS
    "Bjorn Augestad, Meta Systems AS (Metaoptions-0.0.4)",
#endif

#ifdef ABRADFORD
    "Anthony Bradford",
#endif

#ifdef FINRECIPES
    "Bernt Arne Oedegaard (Financial Numerical Recipes in C++)",
#endif

#ifdef SPINSKY
    "Seth Pinsky",
#endif

    NULL
  };

  const gchar *documenters[] = {
    "Anthony Bradford",
    "<info@opensourcefinancialmodels.com>", 
    PACKAGE_URL,
    NULL
  };

  properties->GtkInfo.dialogAbout = gtk_about_dialog_new();

  gtk_window_set_icon(GTK_WINDOW(properties->GtkInfo.dialogAbout), create_pixbuf("images/2.png"));

  gtk_about_dialog_set_program_name(GTK_ABOUT_DIALOG (properties->GtkInfo.dialogAbout), PACKAGE_NAME);

  gtk_about_dialog_set_version(GTK_ABOUT_DIALOG (properties->GtkInfo.dialogAbout), PACKAGE_VERSION);

  gtk_about_dialog_set_comments(GTK_ABOUT_DIALOG (properties->GtkInfo.dialogAbout), 
                                 "The Ultimate Options Calculator\ninfo@opensourcefinancialmodels.com");

  gtk_about_dialog_set_copyright(GTK_ABOUT_DIALOG (properties->GtkInfo.dialogAbout), 
                                  "(C) 2012 Anthony Bradford");

  gtk_about_dialog_set_license(GTK_ABOUT_DIALOG (properties->GtkInfo.dialogAbout), license);
  gtk_about_dialog_set_website(GTK_ABOUT_DIALOG (properties->GtkInfo.dialogAbout), 
                                       PACKAGE_URL);
  gtk_about_dialog_set_website_label(GTK_ABOUT_DIALOG (properties->GtkInfo.dialogAbout), 
                                      "OpenSourceFinancialModels.com");

  gtk_about_dialog_set_authors(GTK_ABOUT_DIALOG(properties->GtkInfo.dialogAbout), authors);
  gtk_about_dialog_set_documenters(GTK_ABOUT_DIALOG(properties->GtkInfo.dialogAbout), documenters);
  gtk_about_dialog_set_translator_credits(GTK_ABOUT_DIALOG(properties->GtkInfo.dialogAbout), "Anthony Bradford\n<info@opensourcefinancialmodels.com>\nhttp://opensourcefinancialmodels.com");

  pngTimer(properties);

  // next line would rotate the images...
  //g_timeout_add(4000, (GSourceFunc) pngTimer, (gpointer) properties);

  gtk_dialog_run(GTK_DIALOG(properties->GtkInfo.dialogAbout));
  gtk_widget_destroy(properties->GtkInfo.dialogAbout);
}

void on_menu_debug_log_activate(GtkWidget *widget, struct _properties *properties)
{
  g_print("on_menu_debug_log_activate()\n");

  if( properties->data.debug == true )
    properties->data.debug = false;
  else
    properties->data.debug = true;
}

void on_menu_text_export_activate(GtkWidget *widget, struct _properties *properties)
{
  g_print("on_menu_text_export_activate()\n");

  properties->textExport = true;
}

void on_menu_source_export_activate(GtkWidget *widget, struct _properties *properties)
{
  g_print("on_menu_source_export_activate()\n");

  // maybe the source should just be pulled from the website
  // with a gtk_show_uri(...) ??

  if( strcmp(option_algorithms[properties->modeltype].sourceCode,"") != 0 )
  {
    g_print("sourceCode exists sourceCode = %s\n", option_algorithms[properties->modeltype].sourceCode );
    display_source(option_algorithms[properties->modeltype].sourceCode, properties);
  }

  if( strcmp(option_algorithms[properties->modeltype].sourceCode2,"") != 0 )
  {
    g_print("sourceCode2 exists sourceCode2 = %s\n", option_algorithms[properties->modeltype].sourceCode2 );
    display_source(option_algorithms[properties->modeltype].sourceCode2, properties);
  }
}

// display_source() fails on displaying metaoptions-0.0.4/src/bisection.c source file
// (under Linux but not Windows) due to some Norwegian in the comments of the source
void display_source(const char *name, const struct _properties *properties)
{
  GtkWidget *window, *scrolled_win, *textview;
  GtkTextBuffer *buffer;

  char windowTitle[400] = { 0 };
  sprintf(windowTitle,"%s Source Export %s", PACKAGE_NAME, name);
  
  window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
  gtk_window_set_title (GTK_WINDOW (window), windowTitle);
  gtk_container_set_border_width (GTK_CONTAINER (window), 10);
  gtk_widget_set_size_request (window, 700, 600);

  gtk_window_set_icon(GTK_WINDOW(window), create_pixbuf("images/2.png"));
  
  textview = gtk_text_view_new();

  FILE *file;
  char *bufferSource;
  unsigned long fileLen;

  //g_print("secondard = src/%s-%s/src/%s\n",PACKAGE_TARNAME,VERSION,name);

  if( access( name, F_OK ) != 0 )
  {
    char name2[200]  = { 0 };
    sprintf(name2,"src/%s-%s/src/%s",PACKAGE_TARNAME,VERSION,name);
    g_print("using second path: %s\n", name2);

    file = fopen(name2, "rb");

  } else
  {
    file = fopen(name, "rb");
  }

  if( !file )
  {
    g_print("display_source(): Unable to open file %s or secondary\n", name);
    return;
  }

  fseek(file, 0, SEEK_END);
  fileLen = ftell(file);
  fseek(file, 0, SEEK_SET);

  bufferSource = (char *) malloc(fileLen+1);

  if( !bufferSource )
  {
    g_print("display_source(): Memory error\n");
    fclose(file);
    return;
  }

  if(!fread(bufferSource, fileLen, 1, file))
  {
    g_print("display_source(): fread() returning 0\n");
  }

  fclose(file);

  bufferSource[fileLen] = 0;
  buffer = gtk_text_view_get_buffer (GTK_TEXT_VIEW (textview));
  
  gtk_text_buffer_set_text (buffer, g_locale_to_utf8 (bufferSource, -1,0,0,0), -1);
  
  scrolled_win = gtk_scrolled_window_new (NULL, NULL);
  gtk_container_add (GTK_CONTAINER (scrolled_win), textview);
  gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (scrolled_win),GTK_POLICY_AUTOMATIC, GTK_POLICY_ALWAYS);
  
  gtk_container_add(GTK_CONTAINER (window), scrolled_win);
  gtk_widget_show_all(window);

  free(bufferSource);
}

void checkForSourceCode(struct _properties *properties)
{
  // maybe the source should just be pulled from the website
  // with a gtk_show_uri(...) ??

  char sourceLocation2[400]  = { 0 };

  sprintf(sourceLocation2,"src/%s-%s/src/%s",PACKAGE_TARNAME,VERSION,option_algorithms[properties->modeltype].sourceCode);
  g_print("sourceLocation2 = %s\n",sourceLocation2);

  if( strcmp( option_algorithms[properties->modeltype].sourceCode,"") != 0 )
  {
    if( access( option_algorithms[properties->modeltype].sourceCode, F_OK ) == 0 )
    {
      g_print("sourceCode is %s\n",option_algorithms[properties->modeltype].sourceCode);
      gtk_widget_show(properties->GtkInfo.menu_source_export);

    } else if( access( sourceLocation2, F_OK ) == 0  )
    {
      g_print("sourceCode is %s\n",sourceLocation2);
      gtk_widget_show(properties->GtkInfo.menu_source_export);

    } else
    {
      g_print("Cannot access sourceCode %s or\n%s\n",option_algorithms[properties->modeltype].sourceCode,sourceLocation2);
      gtk_widget_hide(properties->GtkInfo.menu_source_export);
    }

  } else
  {
    g_print("sourceCode not defined\n");
    gtk_widget_hide(properties->GtkInfo.menu_source_export);
  }
}

void on_menu_ListCategory_activate(GtkWidget *widget, struct _properties *properties)
{
  g_print("on_menu_ListCategory_activate():\n");

  const int num_models = (signed) (sizeofoption_algorithms/sizeof(struct option_algorithm));
  properties->listModelsForGroups = new elementListWithGroup[num_models];

  int index = 0;
  for(index = 0; index < num_models; index++)
  {
    sprintf(properties->listModelsForGroups[index].elementName,"%s",option_algorithms[index].des);
    sprintf(properties->listModelsForGroups[index].groupName,"%s",option_algorithms[index].category);
    properties->listModelsForGroups[index].index = index;
  }

  set_up_combobox_with_array_use_groups(properties->GtkInfo.comboboxModel, properties->listModelsForGroups, num_models, properties->TreeToIndex,1);

}

void on_menu_ListAuthor_activate(GtkWidget *widget, struct _properties *properties)
{
  g_print("on_menu_ListAuthor_activate():\n");

  const int num_models = (signed) (sizeofoption_algorithms/sizeof(struct option_algorithm));

  properties->listModelsForGroups = new elementListWithGroup[num_models];

  int index = 0;
  for(index = 0; index < num_models; index++)
  {
    sprintf(properties->listModelsForGroups[index].elementName,"%s",option_algorithms[index].des);
    sprintf(properties->listModelsForGroups[index].groupName,"%s",option_algorithms[index].source);
    properties->listModelsForGroups[index].index = index;
  }

  set_up_combobox_with_array_use_groups(properties->GtkInfo.comboboxModel, properties->listModelsForGroups, num_models, properties->TreeToIndex,0);

}

void on_menu_List_activate(GtkWidget *widget, struct _properties *properties)
{
  g_print("on_menu_List_activate():\n");

  const int num_models = (signed) (sizeofoption_algorithms/sizeof(struct option_algorithm));
  struct elementList listModels[num_models];

  int index = 0;
  for(index = 0; index < num_models; index++)
  {
    sprintf(listModels[index].elementName,"%d - %s", index + 1, option_algorithms[index].des);
  }

  set_up_combobox_with_array(properties->GtkInfo.comboboxModel, listModels, num_models);
}

void on_menu_settings_activate(GtkWidget *widget, struct _properties *properties)
{
  g_print("on_menu_settings_activate():\n");

  gtk_widget_show_all(properties->GtkInfo.dialogSettings);

  gtk_spin_button_set_value(GTK_SPIN_BUTTON(properties->GtkInfo.spinbuttonSleepDelay),properties->updatedelayseconds);
}
